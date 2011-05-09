
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='http://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
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

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="http://<?php echo $themeroot; ?>images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>

<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<center>
<font size=8 color="darkgreen">caper</font>
<br>
<br>
<font size=5>
<font color="darkgreen">c</font>omparative <font color="darkgreen">a</font>nalyses of <font color="darkgreen">p</font>hylogenetics and <font color="darkgreen">e</font>volution in <font color="darkgreen">R</font>
</font>
</center>
<!-- end of project description -->

<p> The caper package provides a set of tools for conducting phylogenetic comparative analyses in R. The main methods are phylogenetic independent contrasts and phylogenetic generalised least squares methods but the package also provides tools for calculating phylogenetic diversity, examining phylogenetic imbalance and for simulating phylogenies.</p>

<p>
The most recent version of the package can be downloaded from the following <a href="http://r-forge.r-project.org/R/?group_id=983">link</a>. At present, the caper package is only present on R-Forge and is only available for the most recent versions of R. However, the package currently only uses R code and so the following code should install caper on older versions of R:
</p>

<font size=2><tt>install.packages('caper', repos='http://r-forge.r-project.org', type='source', dependencies=TRUE)</tt></font>

<p>The caper package requires the packages 'ape', 'MASS' and 'mvtnorm'. Once caper has been installed, please look at the package vignette for further details:</p>

<font size=2><tt>vignette('caper')</tt></font>

<p>The project development page for the package can be found <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. The package is maintained by David Orme.</p>

</body>
</html>
