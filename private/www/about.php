<?php
include 'include/common.php';
printHeader("Intro to the GATK", "About", "bootstrap_tab_bookmark();");
?>
<script>
    function changeTab( num ) {
    $('#aboutPills li:eq(' + num + ') a').tab('show');
    window.scrollTo(0, document.getElementById('aboutTabs').offsetTop-70);
    }
</script>
<style>
    .space-below {
        margin-bottom : 40px;
    }

    #aboutPills li {
        margin : 0 12px 0 0;
    }

    #aboutPills li a {
        height                : 100px;
        width                 : 100px;
        text-align            : center;
        vertical-align        : middle;
        display               : table-cell;
        padding               : 0;
        border-right          : 100px;
        -moz-border-radius    : 100px !important;
        -webkit-border-radius : 100px !important;
        border-radius         : 100px !important;
    }

    #img-icon {
        height      : 144px;
        width       : 144px;
        float       : left;
        line-height : 144px;
        font-size   : 144px;
    }

    a.nohover:hover {
        text-decoration : none;
    }
</style>

<div class="row">
<nav class="span3">
    <h3>About</h3>
    <ul class="nav nav-list">
        <li class="active"><a href="about">Introduction to the GATK</a></li>
        <li><a href="who-we-are">Who we are</a></li>
        <li><a href="user-stories">User stories</a></li>
    </ul>
</nav>
<!--/.span -->

<div class="span9">

<h1>Introduction to the GATK</h1>
<br/>
<ul id="aboutPills" class="nav nav-pills about">
    <li class="active">
        <a href="#what-is-the-gatk" data-toggle="pill"><h4>What is the<br/>GATK?</h4></a>
    </li>
    <li>
        <a href="#using-the-gatk" data-toggle="pill"><h4>Using the<br/>GATK</h4></a>
    </li>
    <li>
        <a href="#typical-workflows" data-toggle="pill"><h4>Typical<br/>Workflows</h4></a>
    </li>
    <li>
        <a href="#high-performance" data-toggle="pill"><h4>High<br/>Performance</h4></a>
    </li>
    <li>
        <a href="#getting-help" data-toggle="pill"><h4>Getting<br/>Help</h4></a>
    </li>
    <li>
        <a href="#licensing" data-toggle="pill"><h4>Licensing</h4></a>
    </li>
</ul>

<div id="aboutTabs" class="tab-content row-fluid">
<div class="tab-pane fade in active" id="what-is-the-gatk">
    <h2>What is the GATK?<br/>
        <small>Simply what it says on the can: a Toolkit for Genome Analysis</small>
    </h2>
    <hr/>
    <div class="row-fluid">
        <div class="span6">
            <p>Say you have ten exomes and you want to identify the rare mutations they all have in common – the GATK
                can do that. Or you need to know which mutations are specific to a group of patients, as opposed to a
                healthy cohort – the GATK can do that too. In fact, the GATK is the industry standard for such
                analyses.</p>
            <h4>But wait, there’s more!</h4>

            <p>Because of the way it is built, the GATK is highly generic and can be applied to all kinds of datasets
                and genome analysis problems. It can be used for discovery as well as for validation. It’s just as happy
                handling exomes as whole genomes. It can use data generated with a variety of different sequencing
                technologies. And although it was originally developed for human genetics, the GATK has evolved to
                handle genome data from any organism, with any level of ploidy. Your plant has six copies of each
                chromosome? Bring it on.</p>
        </div>
        <div class="span6">
            <div class="thumbnail">
                <img src="img/organisms.png" alt="Organisms">

                <div class="caption">
                    <p>The GATK can handle a variety of organism genomes in addition to humans.</p>
                </div>
            </div>
        </div>
    </div>
    <div class="row-fluid">
        <div class="span6">
            <div class="thumbnail">
                <img src="img/toolduo2.png" alt="Toolbox and toolchain">

                <div class="caption">
                    <p>The toolkit provides a wide set of tools that can be chained into workflows, taking advantage of
                        the common architecture and powerful engine.</p>
                </div>
            </div>
        </div>
        <div class="span6">
            <h4>So what’s in the can?</h4>

            <p>At the heart of the GATK is an industrial-strength infrastructure and engine that handle data access,
                conversion and traversal, as well as high-performance computing features. On top of that lives a rich
                ecosystem of specialized tools, called “walkers”, that you can use out of the box, individually or
                chained into scripted workflows, to perform anything from simple data diagnostics to complex
                “reads-to-results” analyses.</p>

            <p>Some typical workflows are detailed on the next page of this section. Please see the <a href="gatkdocs/">Guides</a>
                section for a complete list of walkers and their capabilities.</p>
        </div>
    </div>
    <hr/>
            <span>
            	<a class="pull-right" href="#using-the-gatk" onclick="changeTab(1)">High Performance <i
                    class="icon-arrow-right"></i></a>
            </span>
</div>

<div class="tab-pane fade" id="using-the-gatk">
    <h2>Using the GATK<br/>
        <small>Get started today</small>
    </h2>
    <hr/>
    <div class="row-fluid">
        <div class="span6">
            <h4>Platform and requirements</h4>

            <p>The GATK is designed to run on Linux and other POSIX-compatible platforms. Yes, that includes MacOS X! If
                you are on any of the above, see <a href="download">the download section</a> for downloading and
                installation instructions. If you’re stuck with Windows, you’re not completely out of luck – see <a
                    href="http://gatk.vanillaforums.com/discussion/7/gatk-on-windows">this post</a> for specific
                instructions on using the GATK with Cygwin. If you’re on something else… no, there are no plans to port
                the GATK to Android or iOS in the near future.</p>

            <p>You will need to have Java installed to run the GATK, and some tools additionally require R to generate
                PDF plots. Version requirements and installation instructions for both are included in the
                platform-specific pages linked to above.</p>
        </div>
        <div class="span6">
            <div class="thumbnail">
                <img src="img/unixes.png" alt="Unix-based operating systems">

                <div class="caption">
                    <p>The GATK is designed to run on Linux and other POSIX-compatible platforms. Yes, that includes
                        MacOS X!</p>
                </div>
            </div>
        </div>
    </div>
    <br/>

    <div class="row-fluid">
        <div class="span6">
            <div class="thumbnail">
                <img src="img/gatkrunsonjava.png" alt="The GATK runs on Java">

                <div class="caption">
                    <p><br>The GATK runs on Java, straight from the command-line.</p>
                </div>
            </div>
        </div>
        <div class="span6">
            <h4>Interface</h4>

            <p>Now here’s the kicker: the GATK does not have a graphical user interface. All tools are called via the
                command-line interface. </p>

            <p>If that is not something you are used to, or you have no idea what that even means, don’t worry. It’s
                easier to learn than you might think, and there are many good <a
                    href="http://lifehacker.com/5633909/who-needs-a-mouse-learn-to-use-the-command-line-for-almost-anything">online
                    tutorials</a> that can get help you get comfortable with the command-line environment. Before you
                know it you’ll be writing scripts to chain tools together into workflows... You don’t need to have any
                programming experience to use the GATK, but you might pick some up along the way!</p>
        </div>
    </div>
    <br/>

    <div class="row-fluid">
        <div class="span6">
            <h4>Command structure and tool arguments</h4>

            <p>All the GATK tools are called using the same basic command structure. Here’s a simple example that counts
                the number of sequence reads in a BAM file:</p>
<pre class="prettyprint linenums">java -jar GenomeAnalysisTK.jar \
-T CountReads \
-R example_reference.fasta \
-I example_reads.bam</pre>
            <p>The <code>-jar</code> argument invokes the GATK engine itself, and the <code>-T</code> argument tells it
                which tool you want to run. Arguments like <code>-R</code> for the genome reference and <code>-I</code>
                for the input file are also given to the GATK engine and can be used with all the tools (see complete
                list of <a href="gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html">available arguments for
                    the GATK engine</a>. Most tools also take additional arguments that are specific to their function.
                These are listed for each tool on that tool’s documentation page, all easily accessible through the <a
                    href="gatkdocs/">Documentation Index</a>.</p>
        </div>
        <div class="span6">
            <div class="thumbnail">
                <img src="img/cmd_output.jpg" alt="Output of the example command">

                <div class="caption">
                    <p>The GATK outputs structured command information, status messages and result summaries to the
                        console.</p>
                </div>
            </div>
        </div>
    </div>
    <br/>

    <p>Please see <a href="http://gatk.vanillaforums.com/discussion/8/video-index">this page</a> for more detailed
        tutorials (including videos!) on using the GATK tools.</p>
    <hr/>
		<span>
			<a class="pull-left" href="#what-is-the-gatk" onclick="changeTab(0)"><i class="icon-arrow-left"></i> What is
                the GATK?</a>
			<a class="pull-right" href="#typical-workflows" onclick="changeTab(2)">Typical Workflows <i
                class="icon-arrow-right"></i></a>
		</span>
</div>

<div class="tab-pane fade" id="typical-workflows">
    <h2>Typical Workflows<br/>
        <small>From sequencing reads to actionable results</small>
    </h2>
    <hr/>
    <p>When you're isolating DNA in the lab, you don't treat the work like isolated, disconnected tasks. Every task is a
        step in a well-documented protocol, carefully developed to optimize yield, purity and to ensure reproducibility
        as well as consistency across all samples and experiments.</p>

    <p>We believe working with NGS data should be exactly the same.</p>

    <p>That's why we have developed industry-standard workflows that are optimized to produce the most accurate results
        from your dataset, with the most efficiency in terms of both manual handling and computational cost.</p>
    <br/>

    <div class="row-fluid">
        <div class="span6">
            <h4>NGS data processing</h4>

            <p>Whatever the sequencing technology you're using, you need to process the raw dataset to make it suitable
                for analysis. This <a href="http://gatk.vanillaforums.com/categories/data-processing">data processing
                    workflow</a> guides you through the necessary steps, with detailed explanations of each operation,
                why it is required and what transformations are applied to the data.</p>
        </div>
        <div class="span6">
            <div class="thumbnail">
                <img src="http://placehold.it/400X240" alt="">
            </div>
        </div>
    </div>
    <br/>

    <div class="row-fluid">
        <div class="span6">
            <h4>Variant discovery, genotyping and filtering</h4>

            <p>Finding sequence variation within and between samples is fairly straightforward. Distinguishing what part
                of that variation is real and assigning the right genotypes is a heck of a lot more difficult. These <a
                    href="http://gatk.vanillaforums.com/categories/variant-discovery">variant discovery</a>, <a
                    href="http://gatk.vanillaforums.com/categories/genotyping">genotyping</a> and <a
                    href="http://gatk.vanillaforums.com/categories/filtering">filtering workflows</a> help you choose
                the parameters that are most appropriate for your dataset and guides you through the necessary steps to
                produce a variant callset that you can trust. Various options are available depending on whether you're
                working with whole genomes or exomes and according to the type, number and coverage depth of your
                samples. [image: variant calling + filtering wf]</p>
        </div>
        <div class="span6">
            <div class="thumbnail">
                <img src="http://placehold.it/400X240" alt="">
            </div>
        </div>
    </div>
    <br/>
    <h4>Best practices for calling variants with the GATK</h4><br/>

    <div class="thumbnail">
        <img src="img/GATK_BPP.png" alt="">

        <div class="caption">
            <p>This reads-to-results variant calling workflow lays out the best practices recommended by our group for
                all the steps involved in calling variants with the GATK. It is used in production at the Broad
                Institute on every genome that rolls out of the sequencing facility. Be sure to check out the series of
                <a href="#">video tutorials</a> devoted to this workflow!</p>
        </div>
    </div>
    <br/>

    <p>Other workflows are available <a href="http://gatk.vanillaforums.com/discussions/tagged/workflow">here</a>
        <del> in the Guides section.</del>
    </p>
    <hr/>
    	<span>
    		<a class="pull-left" href="#using-the-gatk" onclick="changeTab(1)"><i class="icon-arrow-left"></i> Using the
                GATK</a>
        	<a class="pull-right" href="#high-performance" onclick="changeTab(3)">High Performance <i
                class="icon-arrow-right"></i></a>
  		</span>
</div>

<div class="tab-pane fade" id="high-performance">
    <h2>High Performance<br/>
        <small>Built for scalability and parallelism</small>
    </h2>
    <hr/>
    <div class="row-fluid">
        <div class="span6">
            <p>The GATK was built from the ground up with performance in mind.</p>
            <h4>Map/Reduce: it's not just for Google anymore</h4>

            <p>Every GATK walker is built using the Map/Reduce framework, which is basically a strategy to speed up
                performance by breaking down large iterative tasks into shorter segements then merging overall
                results.</p>
            <h4>Muli-threading</h4>

            <p>The GATK takes advantage of the latest processors using multi-threading, i. e. run using multiple cores
                on the same machine, sharing the RAM. To enable multi-threading in the GATK, simply add the <code>-nt
                    x</code> argument to your command line, where <code>x</code> is the number of threads, or cores, you
                want to use.</p>
        </div>
        <div class="span6">
            <div class="thumbnail">
                <img src="img/multithread.png" alt="Multi-threading">

                <div class="caption">
                    <p>The GATK does multi-threading.</p>
                </div>
            </div>
        </div>
    </div>
    <br/>

    <div class="row-fluid">
        <div class="span6">
            <div class="thumbnail">
                <img src="http://placehold.it/400X200" alt="Queue and scatter-gather">

                <div class="caption">
                    <p>Queue manages the entire process, called "scatter-gather", of breaking down big jobs into many
                        smaller ones (scatter) then collecting and merging results when they are done (gather).</p>
                </div>
            </div>
        </div>
        <div class="span6">
            <h4>Out on the farm with Queue</h4>

            <p>Queue is a companion program that allows the GATK to take parallelization to the next level: running jobs
                on a high-performance computing cluster, or server farm.</p>
        </div>
    </div>
    <br/>

    <p>For more information on the high performance features of the GATK and Queue, please see [link].</p>
    <hr/>
		<span>
			<a class="pull-left" href="#typical-workflows" onclick="changeTab(2)"><i class="icon-arrow-left"></i>
                Typical Workflows</a>
			<a class="pull-right" href="#getting-help" onclick="changeTab(4)">Getting Help <i
                class="icon-arrow-right"></i></a>
		</span>
</div>

<div class="tab-pane fade" id="getting-help">
    <h2>Getting Help<br/>
        <small>Fear not! Help is at hand.</small>
    </h2>
    <hr/>
    <p>The GATK has a reputation for being wicked complicated, and it’s not entirely undeserved. With great power comes
        great
        <del>responsibility</del>
        complexity… But fear not! Help is at hand.
    </p>
    <p>[guides (GATK Documentation or G-Docs, tutorials incl. videos) + forum, all integrated/interconnected, all
        searchable by functional category, by keyword, or by browsing the Documentation Index.]</p>
    <br/>
    <h4>GATKDocs</h4>
    <hr/>
    <div class="row-fluid">
        <a class="nohover" href="gatkdocs/"><i id="img-icon" class="icon-book"></i></a>

        <p>This the ultimate reference guide to the GATK. If you’ve ever wondered “Is there a tool that can do X?” or
            “Does this tool have an argument to do Y?”, the <a href="gatkdocs/">GATKDocs</a> should be your first stop.
        </p>

        <p>Every tool in the GATK has its own GATKDoc article detailing what it does and how it does it, as well as all
            the available options, default parameter values and argument names. There are also GATKDocs detailing
            options and arguments of the GATK engine, which are common to all GATK tools, and documenting companion
            software such as Queue.</p>

        <p>It is a living textbook: the articles are generated directly from the source code documentation and updated
            automatically with every new version of the GATK, so you know the information in there is always up to
            date.</p>
    </div>
    <br/>
    <h4>Tutorials and videos</h4>
    <hr/>
    <div class="row-fluid">
        <a class="nohover" href="<?php echo $forum . "discussions/tagged/tutorials" ?>"><i id="img-icon"
                                                                                           class="icon-play-circle"></i></a>

        <p>In addition to the GATKDocs, we are cultivating a growing database of <a
            href="<?php echo $forum . "discussions/tagged/tutorials" ?>">tutorials</a> – some in <a
            href="<?php echo $forum . "discussions/tagged/video" ?>">video</a> form – that will guide you step by step
            through various tasks such as running analyses and troubleshooting errors.</p>

        <p>And because we know there are few things more frustrating than trawling through tutorials that are either too
            basic, too advanced or otherwise not appropriate for your needs, each tutorial is clearly labeled to
            identify the intended audience type ( <span class="label">Analyst</span> or <span
                class="label label-inverse">Developer</span> ), level ( <span class="label label-success">Basic</span>,
            <span class="label label-warning">Intermediate</span>, or <span
                class="label label-important">Advanced</span> ) and prerequisite knowledge.</p>
    </div>
    <br/>
    <h4>Community forum</h4>
    <hr/>
    <div class="row-fluid">
        <a class="nohover" href="<?php echo $forum . "post/discussion" ?>"><i id="img-icon"
                                                                              class="icon-comments-alt"></i></a>

        <p>Finally, if you’ve exhausted all these avenues and still haven't found the answer to your question, </p>

        <h2><a href="<?php echo $forum . "post/discussion" ?>">Ask the forum</a></h2>

        <p>[ask the forum -> user community + GATK team].</p>
    </div>
    <hr/>
        <span>
        	<a class="pull-left" href="#high-performance" onclick="changeTab(3)"><i class="icon-arrow-left"></i> High
                Performance</a>
            <a class="pull-right" href="#licensing" onclick="changeTab(5)">Licensing <i
                class="icon-arrow-right"></i></a>
   		</span>
</div>

<div class="tab-pane fade" id="licensing">
    <h2>Licensing </h2>
    <hr/>
    <div class="row-fluid">
        <div class="span4">
            <h4>Academics</h4>

            <p>Suspendisse vel nunc in purus volutpat dignissim. Nunc id ligula est. Nulla a nisl vitae velit rhoncus
                aliquet in sit amet quam. Ut mollis mollis enim, in mollis orci mollis vel. Ut vel est nec nisl
                fermentum posuere. Vestibulum tellus libero, elementum ac iaculis et, pulvinar non arcu. Maecenas erat
                libero, euismod nec posuere nec, porta ac nunc. Phasellus faucibus, ipsum eu rutrum ultrices, magna mi
                ultricies sapien, ac vestibulum tellus augue non dolor. Phasellus adipiscing, magna eu iaculis blandit,
                augue tellus sagittis magna, sed convallis erat velit id quam. Nam id magna eros.</p>
        </div>
        <div class="span4">
            <h4>Companies</h4>

            <p>Suspendisse vel nunc in purus volutpat dignissim. Nunc id ligula est. Nulla a nisl vitae velit rhoncus
                aliquet in sit amet quam. Ut mollis mollis enim, in mollis orci mollis vel. Ut vel est nec nisl
                fermentum posuere. Vestibulum tellus libero, elementum ac iaculis et, pulvinar non arcu. Maecenas erat
                libero, euismod nec posuere nec, porta ac nunc. Phasellus faucibus, ipsum eu rutrum ultrices, magna mi
                ultricies sapien, ac vestibulum tellus augue non dolor. Phasellus adipiscing, magna eu iaculis blandit,
                augue tellus sagittis magna, sed convallis erat velit id quam. Nam id magna eros.</p>
        </div>
        <div class="span4">
            <h4>Others</h4>

            <p>Suspendisse vel nunc in purus volutpat dignissim. Nunc id ligula est. Nulla a nisl vitae velit rhoncus
                aliquet in sit amet quam. Ut mollis mollis enim, in mollis orci mollis vel. Ut vel est nec nisl
                fermentum posuere. Vestibulum tellus libero, elementum ac iaculis et, pulvinar non arcu. Maecenas erat
                libero, euismod nec posuere nec, porta ac nunc. Phasellus faucibus, ipsum eu rutrum ultrices, magna mi
                ultricies sapien, ac vestibulum tellus augue non dolor. Phasellus adipiscing, magna eu iaculis blandit,
                augue tellus sagittis magna, sed convallis erat velit id quam. Nam id magna eros.</p>
        </div>
    </div>
    <hr/>
      	<span>
      		<a class="pull-left" href="#getting-help" onclick="changeTab(4)"><i class="icon-arrow-left"></i> Getting
                  Help</a>
            <a class="pull-right" href="download">Downloads <i class="icon-arrow-right"></i></a>
   		</span>
</div>
</div>
</div>

</div><!--/row-->

<script>
    function bootstrap_tab_bookmark (selector) { if (selector == undefined) {
    selector = "#aboutPills "; }

    /* Automagically jump on good tab based on anchor */
    $(document).ready(function() {
    url = document.location.href.split('#');
    if(url[1] != undefined) {
    $(selector + '[href=#'+url[1]+']').tab('show');
    }
    });

    var update_location = function (event) {
    document.location.hash = this.getAttribute("href");
    }

    /* Update hash based on tab */
    $(selector + "[data-toggle=pill]").click(update_location);
    $(selector + "[data-toggle=tab]").click(update_location);
    }
</script>

<?php printFooter(); ?>