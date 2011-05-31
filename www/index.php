
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->
<!--
<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
-->
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo 
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>
-->

<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<table>
<tr>
<td><img src=caper64.png></td>
<td><font size=8 color="darkgreen">caper</font></td>
</tr>
<tr>
<td colspan=2>
<font size=5>
<font color="darkgreen" size=5><b>c</b></font>omparative <font color="darkgreen" size=5><b>a</b></font>nalyses of <font color="darkgreen" size=5><b>p</b></font>hylogenetics and <font color="darkgreen" size=5><b>e</b></font>volution in <font color="darkgreen" size=5><b>R</b></font>
</font>
</td>
</tr>
</table>

<br>
<!-- end of project description -->

<p> The caper package provides a set of tools for conducting phylogenetic comparative analyses in R. The main methods are phylogenetic independent contrasts and phylogenetic generalised least squares methods but the package also provides tools for calculating phylogenetic diversity, examining phylogenetic imbalance and for simulating phylogenies.</p>

<h3>Installing <i>caper</i></h3>

<h4>Official versions</h4>
<p>Official 'released' versions of <i>caper</i> for Windows, Mac OS X and Linux are available from the Comprehensive R Archive Network (<a href="http://cran.r.project.org">http://cran.r.project.org</a>) and can be installed from within R using the following code:</p>

<blockquote><font size=2><tt>install.packages('caper')</tt></font></blockquote>

<h4>Development versions</h4>
<p>The package is developed and maintained on the R-Forge website, from which the most recent version of the package can be downloaded or installed. Note that this will include all recent changes and <i>may contain code still in development</i>. Binary packages for Windows and Mac OS X and the most recent version of the source code can be found <a href="http://r-forge.r-project.org/R/?group_id=983">here</a>.</p>

<p>Note that R-Forge only maintains binary packages for the most recent release of R. If you are using the most recent version, then the following code should install the most recent version of <em>caper</em>:</p>

<blockquote><font size=2><tt>install.packages('caper', repos='http://r-forge.r-project.org')</tt></font></blockquote>

<p>If you are using an older version of R, then a binary version of <em>caper</em> is unlikely to be available. Because <em>caper</em> currently only uses pure R code (and not any Fortan or C), it may be possible to install <em>caper</em> on older versions of R using the following code:</p>

<blockquote><font size=2><tt>install.packages('caper', repos='http://r-forge.r-project.org', type='source', dependencies=TRUE)</tt></font></blockquote>

<p>The caper package requires the packages 'ape', 'MASS' and 'mvtnorm'.</p>

<h3>Using <em>caper</em></h3>

<p>In addition to the standard R documentation pages for functions and datasets, the <em>caper</em> package contains two vignettes. The first, the 'caper' vignette, contains a short introduction to phylogenetic comparative methods and step by step guidance on the use of the various functions:</p>

<blockquote><font size=2><tt>vignette('caper')</tt></font></blockquote>

<p>The second vignette contains a set of comparisons of analyses using caper to the same analyses in other software implementations:</p>

<blockquote><font size=2><tt>vignette('caper-benchmarks')</tt></font></blockquote>

<p> Bugs and support requests can be submitted through the <a href="https://r-forge.r-project.org/tracker/?group_id=983">project tracker pages</a>. The project development page for the package can be found <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/">here</a>. The package is maintained by David Orme.</p>

</body>
</html>
