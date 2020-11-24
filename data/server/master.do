cd "D:\Data\workdata\705109"
clear all
set more off

******************************************
** BUILD WORKDATA FROM ADMIN DATA FILES **
******************************************

do rawbuild.do
do deathbuild.do

******************************************
******   MAIN ANALYSIS:             ******
******       1. Empirical results   ******
******       2. Descriptives        ******
******************************************

do main_analysis.do
do descriptives.do

******************************************
******   EXPORT OUTPUT:             ******
******       1. Tables              ******
******       2. Graphs              ******
******       3. Moments             ******
******************************************

do tables.do
do graphs.do
do moments.do
