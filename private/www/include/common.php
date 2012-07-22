<?php
/* Debug printing 
ini_set('display_errors', 'On');
error_reporting(E_ALL | E_STRICT );
*/
$site = "http://gsa-stage:8080/";
$forum = "http://gatk.vanillaforums.com/";
$api = "https://gatk.vanillaforums.com/";

$version = file_get_contents($site . "/gatkdocs/current.version.txt");

if ($version == "") {
    $version = "NOT AVAILABLE";
}

include 'RestReader.php';

function printNav($activeSection)
{
    $site = $GLOBALS['site'];
    $forum = $GLOBALS['forum'];

    echo "<ul class='nav'>";

    // This is the check to get the right nav button on top activated
    echo     "<li", ($activeSection == "Home") ? " class='active'" : "", "><a href='${site}'><i class='icon-home'></i> Home</a></li>";
    echo     "<li", ($activeSection == "About") ? " class='active'" : "", "><a href='${site}about'><i class='icon-info-sign'></i> About</a></li>";
    echo     "<li", ($activeSection == "Guide") ? " class='active'" : "", "><a href='${site}gatkdocs/'><i class='icon-book'></i> Guide</a></li>";
    echo     "<li", ($activeSection == "Community") ? " class='active'" : "", "><a href='${forum}'><i class='icon-comments-alt'></i> Community</a></li>";
    echo     "<li", ($activeSection == "Downloads") ? " class='active'" : "", "><a href='${site}download'><i class='icon-download-alt'></i> Downloads</a></li>";


    echo "</ul>";
}

function printHeader($title = "GATK", $activeSection = "", $onLoad = "")
{
    $site = $GLOBALS['site'];
    ?>

<!DOCTYPE html>
<html lang='en'>
<head>
    <meta charset='utf-8'>
    <title>$title</title>
    <meta name='viewport' content='width=device-width, initial-scale=1.0'>
    <meta name='description' content=''>
    <meta name='author' content=''>

    <!-- Le styles -->

    <link type='text/css' rel='stylesheet' href='${site}css/bootstrap.css'>
    <link type='text/css' rel='stylesheet' href='${site}css/font-awesome.css'>
    <link type='text/css' rel='stylesheet' href='${site}css/prettify.css'>

    <style type='text/css'>
        @media (min-width: 980px) {
            body {
                padding-top    : 80px;
                padding-bottom : 0px;
            }

            .brand {
                padding-top    : 0;
                padding-bottom : 0;
            }

            .nav-collapse {
                margin-top : 18px;
            }
        }

        footer {
            margin-top       : 50px;
            padding-bottom   : 10px;
            padding-top      : 10px;
            background-color : #333;
            color            : #999;
        }
    </style>

    <!-- Le HTML5 shim, for IE6-8 support of HTML5 elements -->
    <!--[if lt IE 9]>
    <script src='http://html5shim.googlecode.com/svn/trunk/html5.js'></script>
    <link rel='stylesheet' href='${site}css/font-awesome-ie7.css'>
    <![endif]-->

    <!-- Le fav and touch icons -->
    <link rel='shortcut icon' href='${site}ico/favicon.ico'>
    <link rel='apple-touch-icon-precomposed' sizes='144x144' href='${site}ico/apple-touch-icon-144-precomposed.png'>
    <link rel='apple-touch-icon-precomposed' sizes='114x114' href='${site}ico/apple-touch-icon-114-precomposed.png'>
    <link rel='apple-touch-icon-precomposed' sizes='72x72' href='${site}ico/apple-touch-icon-72-precomposed.png'>
    <link rel='apple-touch-icon-precomposed' href='${site}ico/apple-touch-icon-57-precomposed.png'>
</head>";
    <?php
    echo "<body onload='prettyPrint();" . $onLoad . "'>";
    echo "<header class='navbar navbar-fixed-top'>
		<div class='navbar-inner'>
			<div class='container'>
				<a class='btn btn-navbar' data-toggle='collapse' data-target='.nav-collapse'>
					<span class='icon-bar'></span>
					<span class='icon-bar'></span>
					<span class='icon-bar'></span>
				</a>
				<a class='brand' href='${site}'><img src='${site}img/gatk-logo-333.png' alt='GATK' style='height:40px;width:93px;' /></a>
				<div class='nav-collapse'>";
    printNav($activeSection);
    echo "<form id='searchform' class='navbar-search pull-right' method='get' action='.$site.'>
                        <label for='s' class='assistive-text hidden'>Search</label>
						<input type='search' class='search-query' name='s' id='s' placeholder='Search' />
					</form>	
				</div><!--/.nav-collapse --></div></div></header>
				<div class='container'>";

}

function printFooter()
{
    $site = $GLOBALS['site'];

    echo "
</div><!--/container-->
<footer>
<div class='container'>
<img src='${site}img/broad-logo-333.png' width='118' height='35' href='http://www.broadinstitute.org/' />
<p>&copy; Broad Institute 2012</p>
</div>
</footer>
<script type='text/javascript' src='${site}js/jquery.js'></script>
<script type='text/javascript' src='${site}js/prettify.js'></script>
<script type='text/javascript' src='${site}js/bootstrap.js'></script>

</body>
</html>";
}

function getForumPosts($tag)
{
    try {
        $restReader = new RestReader($GLOBALS['api'] . "discussions.json/tagged/${tag}");

        $restReader->printAnnouncements(3, "No announcements found with tag ${tag}");
        $restReader->printDiscussions(5, "No posts found with tag ${tag}");
    } catch (Exception $e) {
        printErrorBox($e->getMessage() . " Is the forum up?.");
    }
}

function printErrorBox($message)
{
    echo "<div class='alert alert-error'>
        <button type='button' class='close' data-dismiss='alert'>×</button>
        <strong>Something went wrong!</strong> " . $message . "</div>";
}

?>
