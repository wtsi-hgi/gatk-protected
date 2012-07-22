<?php
include 'include/common.php';
printHeader("Download the GATK", "Downloads", '$("#noAuth").modal({backdrop:"static",keyboard:false}).modal("show");');

if ($_GET['auth'] != "yes") {
    ?>
<style>
    html, body {
        height : 100%;
    }

    .container {
        min-height    : 100%;
        margin-bottom : -330px;
        position      : relative;
    }

    .clearfooter {
        height : 330px;
        clear  : both;
    }

    footer {
        height   : 330px;
        position : relative;
    }
</style>
<div class='modal hide fade' id='noAuth'>
    <div class='modal-header'>
        <button type="button" class="close" data-dismiss="modal">×</button>
        <h3>You must be logged in to download the GATK</h3>
    </div>
    <div class='modal-body'>
        <p>As a replacmeent to our key system, users must register on our forum to download the GATK.</p>
    </div>
    <div class='modal-footer'>
        <a href='#' class='btn'>Register</a>
        <a href='#' class='btn btn-primary'>Login Here &raquo;</a>
    </div>
</div>

<div class="clearfooter"></div>
<script>
    $('#myModal').modal({
    backdrop: 'static',
    keyboard: false
    })
</script>
<?php
    printFooter();
    exit;
}

?>


<style>
    .hero-unit {
        padding : 20px;
    }
</style>
<h2> The current version is <a><?php echo $version; ?></a></h2>
<hr>
<div class="hero-unit">
    <h2>Download GATK 2.0</h2>

    <p>This is the light version of GATK blah blah blah. Look at me while I keep typing Donec id elit non mi porta
        gravida at eget metus. Fusce dapibus, tellus ac cursus commodo, tortor mauris condimentum nibh, ut fermentum
        massa justo sit amet risus. </p>

    <p><a href="auth" class="btn btn-primary btn-large">Proceed <i class="icon-arrow-right icon-white"></i></a></p>
</div>

<div class="hero-unit">
    <h2>Download GATK Lite</h2>

    <p>This is the light version of GATK blah blah blah. Look at me while I keep typing Donec id elit non mi porta
        gravida at eget metus. Fusce dapibus, tellus ac cursus commodo, tortor mauris condimentum nibh, ut fermentum
        massa justo sit amet risus. </p>

    <p><a class="btn btn-large"><i class="icon-github"></i> GitHub &raquo;</a> <a class="btn btn-primary btn-large"><i
        class="icon-download-alt icon-white"></i> Download .jar &raquo;</a></p>
</div>

<div class="hero-unit">
    <p>Download the Bundle </p>

    <p><a class="btn btn-large"><i class="icon-download-alt"></i> Download .tar &raquo;</a></p>
</div>
<?php printFooter(); ?>
