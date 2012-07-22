<?php
include 'include/common.php';
printHeader("GATK | Who we are", "About");
?>
<div class="row">
    <nav class="span3">
        <h3>About</h3>
        <ul class="nav nav-list">
            <li><a href="about">Introduction to the GATK</a></li>
            <li class="active"><a href="who-we-are">Who we are</a></li>
            <li><a href="user-stories">User stories</a></li>
        </ul>
    </nav>
    <!--/.span -->

    <div class="span9">

        <h1>Who we are
            <small>Genome Sequencing and Analysis</small>
        </h1>
        <hr/>
        <p>The <strong>Genome Sequencing and Analysis Group</strong> (GSA) in Medical and Population Genetics at the
            Broad Institute is a team of computational biologists, software engineers, and hosted students and
            researchers developing algorithms for next generation DNA sequencers for medical and population genetics and
            cancer applications, as well as applying these algorithms to answer fundamental scientific questions. GSA
            has extensive experience with processing of next-generation DNA sequencer data as well as genotyping and
            validation data along with downstream analysis of this data for medical and population genetics studies. The
            method development arm of GSA has created a powerful framework in the The Genome Analysis Toolkit for
            analysis of next-generation sequencing data and analysis of variation discovered by NGS. These tools are now
            widely used in many NGS projects, including the 1000 Genomes Project, <a rel="nofollow"
                                                                                     href="http://cancergenome.nih.gov/">The
                Cancer Genome Atlas</a>, the Broad's production sequencing pipeline, as well as at many other sequencing
            centers and individual labs with sequencing machines. </p>

        <div class="row-fluid">
            <div class="span6">
                <h2>Group members</h2>
                <ul>
                    <li>Group Manager
                        <ul>
                            <li>Mark A. DePristo, Associate Director of MPG Analysis</li>
                        </ul>
                    </li>
                </ul>
                <ul>
                    <li>Development Team
                        <ul>
                            <li>Eric Banks, Team Lead</li>
                            <li>Guillermo del Angel</li>
                            <li>Ryan Poplin</li>
                            <li>Chris Hartl</li>
                        </ul>
                    </li>
                </ul>
                <ul>
                    <li>Technology Development Team
                        <ul>
                            <li>Mauricio Carneiro, Team Lead</li>
                            <li>Geraldine Van der Auwera</li>
                            <li>Roger Zurawicki</li>
                        </ul>
                    </li>
                </ul>
                <ul>
                    <li>Production data processing
                        <ul>
                            <li>Khalid Shakir, Team Lead</li>
                            <li>David Roazen</li>
                            <li>Joel Thibault</li>
                        </ul>
                    </li>
                </ul>
            </div>
            <div class="span6">
                <div class="thumbnail">
                    <img src="http://placehold.it/600X480" alt="">

                    <div class="caption">
                        <h5>Thumbnail label</h5>

                        <p>Cras justo odio, dapibus ac facilisis in, egestas eget quam. Donec id elit non mi porta
                            gravida at eget metus. Nullam id dolor id nibh ultricies vehicula ut id elit.</p>
                    </div>
                </div>
            </div>
        </div>
    </div>
</div><!--/row-->
<?php printFooter(); ?>