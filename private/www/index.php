<?php
include 'include/common.php';
printHeader("GATK Main Page", "Home");
?>
<style>
    div.big-link a p {
        color : #000;
    }

    img.visible-desktop {
        padding : 5px;
    }

    .version {
        font-size : large;
    }

    div.big-link a:hover {
        text-decoration : none;
    }

    .btn-flat {
        background : #269ABC !important;
    }

    .btn-flat:active {
        background : transparent;
    }
</style>
<div class="hero-unit">
    <img src="http://placehold.it/250x200" align="right">

    <p>The GATK is a structured software library that makes writing efficient analysis tools using next-generation
        sequencing data very easy, and second it's a suite of tools for working with human medical resequencing projects
        such as 1000 Genomes and The Cancer Genome Atlas. These tools include things like a depth of coverage analyzers,
        a quality score recalibrator, a SNP/indel caller and a local realigner.</p>

    <p><a href="about" class="btn btn-primary btn-large btn-flat">Learn more &raquo;</a></p>
</div>
<br/>
<div class="row">
    <div class="span8">
        <div class="row-fluid big-link">
            <a href="about">
                <div class="span6"><img class="visible-desktop" src="img/circle-icon.png" alt="" align="left"/>

                    <h3>About</h3>

                    <p>Quick description</p></div>
            </a>
            <a href="gatkdocs/">
                <div class="span6"><img class="visible-desktop" src="img/circle-icon.png" alt="" align="left"/>

                    <h3>Guides</h3>

                    <p>Quick description</p></div>
            </a>
        </div>
        <br/>

        <div class="row-fluid big-link">
            <a href="vanilla/">
                <div class="span6"><img class="visible-desktop" src="img/circle-icon.png" alt="" align="left"/>

                    <h3>Community</h3>

                    <p>Quick description</p></div>
            </a>
            <a href="http://gatk.vanillaforums.com/discussion/8/video-index">
                <div class="span6"><img class="visible-desktop" src="img/circle-icon.png" alt="" align="left"/>

                    <h3>Videos</h3>

                    <p>Quick description</p></div>
            </a>
        </div>
    </div>

    <div class="span4">
        <div class="well">
            <h3>Latest stable release</h3>
            <span class="lead"><?php echo $version; ?></span>

            <p><a href="#">Release Notes</a> (2012-07-11)</p>

            <p style="text-align: center;"><a href="download" class="btn btn-primary btn-flat">Download now <i
                class="icon-external-link"></i></a></p>
        </div>
    </div>
    <!--/span-->
</div>

<br/>

</div>

<?php printFooter(); ?>