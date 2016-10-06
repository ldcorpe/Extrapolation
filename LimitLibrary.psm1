# Functions to help with limit setting

function Invoke-LumiRun ($inputFile, $abcdInfo, $lum, $abcdError, $jobID, $dataset) {
	$stubname = "LimitFinderLog-$lum-$abcdError-$jobID-$dataset-$abcdInfo".Replace(" ", "-")
	$logfile = "$stubname.log".Replace(" ", "-")
	if (-not $(Test-Path $logfile)) {
		$r = .\Release\ExtrapLimitFinder.exe -extrapFile $inputFile.FullName -UseAsym -nA $abcdInfo[0] -nB $abcdInfo[1] -nC $abcdInfo[2] -nD $abcdInfo[3] -OutputFile "$stubname.root" -Luminosity $lum -ABCDError $abcdError 2>&1
		$r | Set-Content $logfile
	}

	$limitExp = Get-TH1Property "$stubname.root" Minimum xsec_BR_95CL
	$limitObs = Get-TH1Property "$stubname.root" Minimum xsec_BR_events__limit
	return ($limitExp, $limitObs)
}

# What we want others to see
Export-ModuleMember -Function Invoke-LumiRun
