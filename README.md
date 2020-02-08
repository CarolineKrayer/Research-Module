## Project workspace Research Module in Econometrics and Statistics

This repository forms the workspace of our project for the Research Module in
Econometrics and Statistics at the University of Bonn (Winter Semester 2019/2020).
The project is written entirely in R.

*Title*: Regression Discontinuity Design: A Test for Manipulation of the Running Variable

*Instructor*: JProf. Dr. Dominik Liebl

*Authors*: Sofia Badini, Max Schäfer, Caroline Krayer


<hr />



#### Abstract
In this paper we deal with the manipulation test by McCrary (2008) which is commonly used in Regression Discontinuity Design (RDD) applications to detect manipulation of the running variable. To investigate the test's ability to reject the null-hypothesis and highlight the test's properties, we design a simulation study where we draw a running variable from various distributions and apply a selection rule that creates a discontinuity in the underlying density function. In line with derived asymptotic properties, the test's power increases both with the number of observations and the gap size. Further, we investigate the quality of the standard normal limiting distribution for different sample sizes and different degrees of undersmoothing. We show that normality improves with an increasing sample size but the influence of undersmoothing is inconclusive. However, undersmoothing the bandwidth improves the centering of the test statistic's distribution in general, but does not necessarily result in the correct centering at zero. Our results suggest that the test's size is robust to bandwidth choice. Moreover, to study manipulation of social program eligibility in Colombia we apply the test to a large administrative dataset where the running variable is already discrete. We find that the test rejects for all years in which we strongly believe that manipulation is present, while for years with no further evidence of manipulation the test seems to reject the null-hypothesis too often. In the end, both the simulation study and the application to real data highlight that parameter choice is non-trivial and the researcher should not rely on rules-of-thumb alone. The test's performance depends on the data structure at hand and further investigation of manipulation is always required. Despite these practical limitations, our simulation study and application to real data showcase the relevance for this test as an important source of information when evaluating the validity of a RDD.


<hr />


We present the paper, the code and the 25-min presentation held in class. Please feel free to ask any question, or to comment on the work.



[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE.txt)
