package org.broadinstitute.sting.queue.engine

import org.jgrapht.traverse.TopologicalOrderIterator
import org.jgrapht.graph.SimpleDirectedGraph
import scala.collection.JavaConversions
import scala.collection.JavaConversions._
import org.jgrapht.alg.CycleDetector
import org.jgrapht.EdgeFactory
import org.jgrapht.ext.DOTExporter
import java.io.File
import org.jgrapht.event.{TraversalListenerAdapter, EdgeTraversalEvent}
import org.broadinstitute.sting.queue.{QSettings, QException}
import org.broadinstitute.sting.queue.function.scattergather.{GatherFunction, ScatterGatherableFunction}
import org.broadinstitute.sting.queue.function.{InProcessFunction, CommandLineFunction, QFunction}
import org.broadinstitute.sting.queue.util.{JobExitException, LsfKillJob, Logging}

/**
 * The internal dependency tracker between sets of function input and output files.
 */
class QGraph extends Logging {
  var dryRun = true
  var bsubAllJobs = false
  var skipUpToDateJobs = false
  var dotFile: File = _
  var expandedDotFile: File = _
  var qSettings: QSettings = _
  var debugMode = false
  private val jobGraph = newGraph

  /**
   * Adds a QScript created CommandLineFunction to the graph.
   * @param command Function to add to the graph.
   */
  def add(command: QFunction) {
    try {
      command.qSettings = this.qSettings
      command.freeze
      addEdge(new FunctionEdge(command))
    } catch {
      case e: Exception =>
        throw new QException("Error adding function: " + command, e)
    }
  }


  private def scatterGatherable(edge: FunctionEdge) = {
    edge.function match {
      case scatterGather: ScatterGatherableFunction if (scatterGather.scatterGatherable) => true
      case _ => false
    }
  }


  /**
   * Checks the functions for missing values and the graph for cyclic dependencies and then runs the functions in the graph.
   */
  def run = {
    val numMissingValues = fillGraph
    val isReady = numMissingValues == 0

    if (isReady || this.dryRun) {
      runJobs()
    }

    if (numMissingValues > 0) {
      logger.error("Total missing values: " + numMissingValues)
    }

    if (isReady && this.dryRun) {
      logger.info("Dry run completed successfully!")
      logger.info("Re-run with \"-run\" to execute the functions.")
    }
  }

  private def fillGraph = {
    fill
    if (dotFile != null)
      renderToDot(dotFile)
    var numMissingValues = validate

    if (numMissingValues == 0 && bsubAllJobs) {
      logger.debug("Scatter gathering jobs.")
      var scatterGathers = List.empty[FunctionEdge]
      loop({
        case edge: FunctionEdge if (scatterGatherable(edge)) =>
          scatterGathers :+= edge
      })

      var addedFunctions = List.empty[QFunction]
      for (scatterGather <- scatterGathers) {
        val functions = scatterGather.function.asInstanceOf[ScatterGatherableFunction].generateFunctions()
        if (this.debugMode)
          logger.debug("Scattered into %d parts: %n%s".format(functions.size, functions.mkString("%n".format())))
        addedFunctions ++= functions
      }

      this.jobGraph.removeAllEdges(scatterGathers)
      prune
      addedFunctions.foreach(this.add(_))

      fill
      val scatterGatherDotFile = if (expandedDotFile != null) expandedDotFile else dotFile
      if (scatterGatherDotFile != null)
        renderToDot(scatterGatherDotFile)
      numMissingValues = validate
    }

    numMissingValues
  }

  def checkStatus = {
    // build up the full DAG with scatter-gather jobs
    fillGraph
    logStatus
  }

  /**
   * Walks up the graph looking for the previous command line edges.
   * @param function Function to examine for a previous command line job.
   * @param qGraph The graph that contains the jobs.
   * @return A list of prior jobs.
   */
  def previousFunctions(edge: QEdge) : List[FunctionEdge] = {
    var previous = List.empty[FunctionEdge]

    val source = this.jobGraph.getEdgeSource(edge)
    for (incomingEdge <- this.jobGraph.incomingEdgesOf(source)) {
      incomingEdge match {

      // Stop recursing when we find a job along the edge and return its job id
        case functionEdge: FunctionEdge => previous :+= functionEdge

        // For any other type of edge find the jobs preceding the edge
        case edge: QEdge => previous ++= previousFunctions(edge)
      }
    }
    previous
  }

  /**
   * Fills in the graph using mapping functions, then removes out of date
   * jobs, then cleans up mapping functions and nodes that aren't need.
   */
  private def fill = {
    fillIn
    prune
  }

  /**
   * Looks through functions with multiple inputs and outputs and adds mapping functions for single inputs and outputs.
   */
  private def fillIn = {
    // clone since edgeSet is backed by the graph
    JavaConversions.asSet(jobGraph.edgeSet).clone.foreach {
      case cmd: FunctionEdge => {
        addCollectionOutputs(cmd.outputs)
        addCollectionInputs(cmd.inputs)
      }
      case map: MappingEdge => /* do nothing for mapping edges */
    }
  }

  private def getReadyJobs = {
    var readyJobs = List.empty[FunctionEdge]
    loop({
      case f: FunctionEdge => {
        if (this.previousFunctions(f).forall(_.status == RunnerStatus.DONE) && f.status == RunnerStatus.PENDING)
          readyJobs :+= f
      }
    })
    readyJobs
  }

  private def getRunningJobs = {
    var runningJobs = List.empty[FunctionEdge]
    loop({
      case f: FunctionEdge => {
        if (f.status == RunnerStatus.RUNNING)
          runningJobs :+= f
      }
    })
    runningJobs
  }

  /**
   *  Removes mapping edges that aren't being used, and nodes that don't belong to anything.
   */
  private def prune = {
    var pruning = true
    while (pruning) {
      pruning = false
      val filler = jobGraph.edgeSet.filter(isFiller(_))
      if (filler.size > 0) {
        jobGraph.removeAllEdges(filler)
        pruning = true
      }
    }

    jobGraph.removeAllVertices(jobGraph.vertexSet.filter(isOrphan(_)))
  }

  /**
   * Validates that the functions in the graph have no missing values and that there are no cycles.
   * @return Number of missing values.
   */
  private def validate = {
    var numMissingValues = 0
    JavaConversions.asSet(jobGraph.edgeSet).foreach {
      case cmd: FunctionEdge =>
        val missingFieldValues = cmd.function.missingFields
        if (missingFieldValues.size > 0) {
          numMissingValues += missingFieldValues.size
          logger.error("Missing %s values for function: %s".format(missingFieldValues.size, cmd.function.description))
          for (missing <- missingFieldValues)
            logger.error("  " + missing)
        }
      case map: MappingEdge => /* do nothing for mapping edges */
    }

    val detector = new CycleDetector(jobGraph)
    if (detector.detectCycles) {
      logger.error("Cycles were detected in the graph:")
      for (cycle <- detector.findCycles)
        logger.error("  " + cycle)
      throw new QException("Cycles were detected in the graph.")
    }

    numMissingValues
  }

  /**
   * Runs the jobs by traversing the graph.
   */
  private def runJobs() = {
    loop({ case f: FunctionEdge => {
      val isDone = this.skipUpToDateJobs &&
              f.status == RunnerStatus.DONE &&
              this.previousFunctions(f).forall(_.status == RunnerStatus.DONE)
      if (!isDone)
        f.resetPending()
    }})

    var readyJobs = getReadyJobs
    var runningJobs = Set.empty[FunctionEdge]
    while (readyJobs.size + runningJobs.size > 0) {
      var exitedJobs = List.empty[FunctionEdge]
      runningJobs.foreach(runner => {
        if (runner.status != RunnerStatus.RUNNING)
          exitedJobs :+= runner
      })
      exitedJobs.foreach(runner => runningJobs -= runner)

      readyJobs.foreach(f => {
        f.runner = newRunner(f.function)
        f.runner.start()
        if (f.status == RunnerStatus.RUNNING) {
          runningJobs += f
        }
      })

      if (readyJobs.size == 0 && runningJobs.size > 0)
        Thread.sleep(30000L)
      readyJobs = getReadyJobs
    }
  }

  private def newRunner(f: QFunction) = {
    if (this.dryRun)
      new DryRunner(f)
    else {
      f match {
        case cmd: CommandLineFunction =>
          if (this.bsubAllJobs)
            new LsfJobRunner(cmd)
          else
            new ShellJobRunner(cmd)
        case inProc: InProcessFunction =>
          new InProcessRunner(inProc)
        case _ =>
          throw new QException("Unexpected function: " + f)
      }
    }
  }

  /**
   * Tracks analysis status.
   */
  private class AnalysisStatus(val analysisName: String) {
    var status = RunnerStatus.PENDING
    var scatter = new ScatterGatherStatus
    var gather = new ScatterGatherStatus
  }

  /**
   * Tracks scatter gather status.
   */
  private class ScatterGatherStatus {
    var total = 0
    var done = 0
    var failed = 0
  }

  /**
   * Gets job statuses by traversing the graph and looking for status-related files
   */
  private def logStatus = {
    var statuses = Map.empty[String, AnalysisStatus]
    loop({
      case edgeCLF: FunctionEdge if (edgeCLF.function.analysisName != null) =>
        updateStatus(statuses.get(edgeCLF.function.analysisName) match {
          case Some(status) => status
          case None =>
            val status = new AnalysisStatus(edgeCLF.function.analysisName)
            statuses += edgeCLF.function.analysisName -> status
            status
        }, edgeCLF)
    })

    statuses.values.toList.sortBy(_.analysisName).foreach(status => {
      if (status.scatter.total + status.gather.total > 0) {
        var sgStatus = RunnerStatus.PENDING
        if (status.scatter.failed + status.gather.failed > 0)
          sgStatus = RunnerStatus.FAILED
        else if (status.scatter.done + status.gather.done == status.scatter.total + status.gather.total)
          sgStatus = RunnerStatus.DONE
        else if (status.scatter.done + status.gather.done > 0)
          sgStatus = RunnerStatus.RUNNING
        status.status = sgStatus
      }

      var info = status.analysisName + ": [" + status.status + "]"
      if (status.scatter.total + status.gather.total > 1) {
        info += formatSGStatus(status.scatter, "s")
        info += formatSGStatus(status.gather, "g")
      }
      logger.info(info)
    })
  }

  /**
   * Updates a status map with scatter/gather status information (e.g. counts)
   */
  private def updateStatus(stats: AnalysisStatus, edge: FunctionEdge) = {
    if (edge.function.isInstanceOf[GatherFunction]) {
      updateSGStatus(stats.gather, edge)
    } else if (edge.function.isInstanceOf[ScatterGatherableFunction]) {
      updateSGStatus(stats.scatter, edge)
    } else {
      stats.status = edge.status
    }
  }

  private def updateSGStatus(stats: ScatterGatherStatus, edge: FunctionEdge) = {
    stats.total += 1
    edge.status match {
      case RunnerStatus.DONE => {
        stats.done += 1
      }
      case RunnerStatus.FAILED => {
        stats.failed += 1
      }
      /* can't tell the difference between pending and running right now! */
      case RunnerStatus.PENDING =>
      case RunnerStatus.RUNNING =>
    }
  }

  /**
   * Formats a status into nice strings
   */
  private def formatSGStatus(stats: ScatterGatherStatus, prefix: String) = {
    " %s:%dt/%dd/%df".format(
      prefix, stats.total, stats.done, stats.failed)
  }

  /**
   *   Creates a new graph where if new edges are needed (for cyclic dependency checking) they can be automatically created using a generic MappingFunction.
   * @return A new graph
   */
  private def newGraph = new SimpleDirectedGraph[QNode, QEdge](new EdgeFactory[QNode, QEdge] {
    def createEdge(input: QNode, output: QNode) = new MappingEdge(input.files, output.files)})

  private def addEdge(edge: QEdge) = {
    val inputs = QNode(edge.inputs)
    val outputs = QNode(edge.outputs)
    val newSource = jobGraph.addVertex(inputs)
    val newTarget = jobGraph.addVertex(outputs)
    val removedEdges = jobGraph.removeAllEdges(inputs, outputs)
    val added = jobGraph.addEdge(inputs, outputs, edge)
    if (this.debugMode) {
      logger.debug("Mapped from:   " + inputs)
      logger.debug("Mapped to:     " + outputs)
      logger.debug("Mapped via:    " + edge)
      logger.debug("Removed edges: " + removedEdges)
      logger.debug("New source?:   " + newSource)
      logger.debug("New target?:   " + newTarget)
      logger.debug("")
    }
  }

  /**
   * Checks to see if the set of files has more than one file and if so adds input mappings between the set and the individual files.
   * @param files Set to check.
   */
  private def addCollectionInputs(files: Set[File]): Unit = {
    if (files.size > 1)
      for (file <- files)
        addMappingEdge(Set(file), files)
  }

  /**
   * Checks to see if the set of files has more than one file and if so adds output mappings between the individual files and the set.
   * @param files Set to check.
   */
  private def addCollectionOutputs(files: Set[File]): Unit = {
    if (files.size > 1)
      for (file <- files)
        addMappingEdge(files, Set(file))
  }

  /**
   * Adds a directed graph edge between the input set and the output set if there isn't a direct relationship between the two nodes already.
   * @param input Input set of files.
   * @param output Output set of files.
   */
  private def addMappingEdge(input: Set[File], output: Set[File]) = {
    val hasEdge = input == output ||
            jobGraph.getEdge(QNode(input), QNode(output)) != null ||
            jobGraph.getEdge(QNode(output), QNode(input)) != null
    if (!hasEdge)
      addEdge(new MappingEdge(input, output))
  }

  /**
   * Returns true if the edge is mapping edge that is not needed because it does
   * not direct input or output from a user generated CommandLineFunction.
   * @param edge Edge to check.
   * @return true if the edge is not needed in the graph.
   */
  private def isFiller(edge: QEdge) = {
    if (edge.isInstanceOf[MappingEdge]) {
      if (jobGraph.outgoingEdgesOf(jobGraph.getEdgeTarget(edge)).size == 0)
        true
      else if (jobGraph.incomingEdgesOf(jobGraph.getEdgeSource(edge)).size == 0)
        true
      else false
    } else false
  }

  /**
   * Returns true if the node is not connected to any edges.
   * @param node Node (set of files) to check.
   * @return true if this set of files is not needed in the graph.
   */
  private def isOrphan(node: QNode) =
    (jobGraph.incomingEdgesOf(node).size + jobGraph.outgoingEdgesOf(node).size) == 0

  /**
   * Utility function for looping over the internal graph and running functions.
   * @param edgeFunction Optional function to run for each edge visited.
   * @param nodeFunction Optional function to run for each node visited.
   */
  private def loop(edgeFunction: PartialFunction[QEdge, Unit] = null, nodeFunction: PartialFunction[QNode, Unit] = null) = {
    val iterator = new TopologicalOrderIterator(this.jobGraph)
    iterator.addTraversalListener(new TraversalListenerAdapter[QNode, QEdge] {
      override def edgeTraversed(event: EdgeTraversalEvent[QNode, QEdge]) = event.getEdge match {
        case cmd: FunctionEdge => if (edgeFunction != null && edgeFunction.isDefinedAt(cmd)) edgeFunction(cmd)
        case map: MappingEdge => /* do nothing for mapping functions */
      }
    })
    iterator.foreach(node => if (nodeFunction != null && nodeFunction.isDefinedAt(node)) nodeFunction(node))
  }

  /**
   * Outputs the graph to a .dot file.
   * http://en.wikipedia.org/wiki/DOT_language
   * @param file Path to output the .dot file.
   */
  private def renderToDot(file: java.io.File) = {
    val out = new java.io.FileWriter(file)

    // todo -- we need a nice way to visualize the key pieces of information about commands.  Perhaps a
    // todo -- visualizeString() command, or something that shows inputs / outputs
    val ve = new org.jgrapht.ext.EdgeNameProvider[QEdge] {
        def getEdgeName(function: QEdge) = if (function.dotString == null) "" else function.dotString.replace("\"", "\\\"")
    }

    //val iterator = new TopologicalOrderIterator(qGraph.jobGraph)
    (new DOTExporter(new org.jgrapht.ext.IntegerNameProvider[QNode](), null, ve)).export(out, jobGraph)

    out.close
  }

  /**
   * Returns true if any of the jobs in the graph have a status of failed.
   * @return true if any of the jobs in the graph have a status of failed.
   */
  def hasFailed = {
    this.jobGraph.edgeSet.exists(edge => {
      edge.isInstanceOf[FunctionEdge] && edge.asInstanceOf[FunctionEdge].status == RunnerStatus.FAILED
    })
  }

  /**
   * Kills any forked jobs still running.
   */
  def shutdown() {
    val lsfJobs = getRunningJobs.filter(_.runner.isInstanceOf[LsfJobRunner]).map(_.runner.asInstanceOf[LsfJobRunner].job)
    if (lsfJobs.size > 0) {
      for (jobs <- lsfJobs.grouped(10)) {
        try {
          val bkill = new LsfKillJob(jobs)
          logger.info(bkill.command)
          bkill.run()
        } catch {
          case jee: JobExitException =>
            logger.error("Unable to kill all jobs:%n%s".format(jee.getMessage))
          case e =>
            logger.error("Unable to kill jobs.", e)
        }
      }
    }
  }
}
