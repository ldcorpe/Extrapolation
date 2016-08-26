#
# Simple script to invoke a double extrapolation test. Results are dumped to root files and to the terminal window.
#

Param(
	[Parameter(Mandatory=$True,Position=1)]
	[string]$dsname,
	[Parameter(Mandatory=$True,Position=2)]
	[string]$ctau5m,
	[Parameter(Mandatory=$True,Position=3)]
	[string]$ctau9m
)

# Just run it all four times. This takes quite some time!

./Release/ExtrapolateByBeta.exe --muonTreeFile "C:\Users\gordo\Downloads\LLPExtrapolationMCTree-$($dsname)lt5mW.root" --output "extrapolate-$($dsname)lt5mw-flat.root" --ctau $ctau5m -UseFlatBeta
Write-Host ""
./Release/ExtrapolateByBeta.exe --muonTreeFile "C:\Users\gordo\Downloads\LLPExtrapolationMCTree-$($dsname)lt9mW.root" --output "extrapolate-$($dsname)lt9mw-flat.root" --ctau $ctau9m -UseFlatBeta
Write-Host ""
Write-Host ""

./Release/ExtrapolateByBeta.exe --muonTreeFile "C:\Users\gordo\Downloads\LLPExtrapolationMCTree-$($dsname)lt5mW.root" --output "extrapolate-$($dsname)lt5mw-beta.root" --ctau $ctau5m
Write-Host ""
./Release/ExtrapolateByBeta.exe --muonTreeFile "C:\Users\gordo\Downloads\LLPExtrapolationMCTree-$($dsname)lt9mW.root" --output "extrapolate-$($dsname)lt9mw-beta.root" --ctau $ctau9m
