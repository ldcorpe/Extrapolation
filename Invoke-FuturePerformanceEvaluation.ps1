#
# This will run the official limit several times with different configurations.
#

#
# Config the default settings
#

# The number of events in each of the A, B, C, and D regions
$abcd = (24, 16, 39, 34)

# Lumi in fb of the current results
$lumi = 3.2
$otherLumis = (30, 300)

# the size of the ABCD error
$abcdError = 0.36
$otherAbcdErrors = (0.10, 0.20)

# Fraction of background we find.
$abcdFactors = (0.90, 0.80, 0.60)

# What samples to look at and where to get the MC root files.
#$singalSamples = (("1000pi400lt5m", 104), ("400pi50lt5m", 116))
$singalSamples = (("1000pi400lt5m", 104), ("400pi50lt5m", 116))

###################
# We need some help to get at some ROOT files... So we are going to have to pull a slippery one here since we aren't quite all the way
# there with how to distribute this.

# Make sure we are runnign the 32 bit of powershell, otherwise we can't
# load in the ROOT versions of things
if ($env:PROCESSOR_ARCHITECTURE -ne "x86") {
	Write-Error "You are not running the 32 bit version of powershell. Do it by starting powershell with 'C:\Windows\SysWOW64\WindowsPowerShell\v1.0\powershell.exe'"
	return
}

# Add ROOT to the path if not there already.
if (-not $($env:PATH.Contains(("5.34.36")))) {
	$env:ROOTSYS="C:\Users\Gordon\AppData\Local\Temp\root\root-5.34.36-vc12"
	$env:PATH="$env:PATH;$env:ROOTSYS\bin"
}

# Get the ROOT access guy
$pathToCMD = $PWD.Path + "\..\..\PSCernRootUtilities\PSCernRootCommands\bin\x86\Debug\PSCernRootCommands.dll"
$fileio = [System.IO.FileInfo] $pathToCMD
if (-Not $($fileio.Exists)) {
	$msg = "Unable to locate PSCernRootCommands at " + $fileio.FullName
	Write-Error $msg
	return
}
Import-Module $fileio.FullName

#PS C:\Users\Gordon> $env:ROOTSYS="C:\Users\Gordon\AppData\Local\Temp\root\root-5.34.36-vc12"
#PS C:\Users\Gordon> $env:PATH="$env:PATH;$env:ROOTSYS\bin"

########################

function Invoke-LumiRun ($inputFile, $abcdInfo, $lum, $abcdError, $jobID, $dataset) {
	$stubname = "LimitFinderLog-$lum-$abcdError-$jobID-$dataset-$abcdInfo".Replace(" ", "-")
	$logfile = "$stubname.log".Replace(" ", "-")
	if (-not $(Test-Path $logfile)) {
		$r = .\Release\ExtrapLimitFinder.exe -extrapFile $inputFile.FullName -UseAsym -nA $abcdInfo[0] -nB $abcdInfo[1] -nC $abcdInfo[2] -nD $abcdInfo[3] -OutputFile "$stubname.root" -Luminosity $lum -ABCDError $abcdError 2>1
		$r | Set-Content $logfile
	}

	$limitExp = Get-TH1Property "$stubname.root" Minimum xsec_BR_95CL
	$limitObs = Get-TH1Property "$stubname.root" Minimum xsec_BR_events__limit
	return ($limitExp, $limitObs)
}

function Get-ABCDScaled ($abcdOriginal, $scale) {
	$a = $abcdOriginal[0] * $scale
	$b = $abcdOriginal[1] * $scale
	$c = $abcdOriginal[2] * $scale
	$d = $abcdOriginal[3] * $scale
	return ($a, $b, $c, $d)
}

# Function that will run a start job and setup the output to be collected.
$jobs = @()
function Start-LumiCalculation ($message, $inputFile, $abcdInfo, $lum, $abcdError, $jobID, $dataset) {
	# Rate limit
	while (@(Get-Job -State Running).Count -ge 3) {
		Start-Sleep -Seconds 10
	}

	$jobs += Start-Job -ScriptBlock {
		Param($message, $inputFile, $abcdInfo, $lum, $abcdError, $jobID, $dataset)
		$limit = Invoke-LumiRun($inputFile, $abcdInfo, $lum, $abcdError, $jobID, $dataset)
		return ($limit, $dataset, $message)
	}
}

#########
# Run it!

Write-Output "Dataset, Extrapolation, Expected, Observed"
foreach ($signal in $singalSamples) {
	$dataset = $signal[0]
	$jobID = $signal[1]

	# Get the input file from jenkins
	$inputFile = Get-JenkinsArtifact -JobUri http://jenks-higgs.phys.washington.edu:8080/view/LLP/job/Limit-LifetimeEfficiencies/ -JobId $jobID -ArtifactName lifetime_eff_$dataset.root

	# Get the default running result
	Start-LumiCalculation "Default" $inputFile $abcd $lumi $abcdError $jobID $dataset
	#Write-Output "$dataset, Default, $($limit[0]), $($limit[1])"

	# Next, do the increased amount of lumi
	foreach ($newLumi in $otherLumis) {
		$factor = $newLumi / $lumi
		$newabcd = Get-ABCDScaled $abcd $factor
	    Start-LumiCalculation "$newLumi fb" $inputFile $newabcd $newLumi $abcdError $jobID $dataset
		#Write-Output "$dataset, $newLumi fb, $($limit[0]), $($limit[1])"
	}

	# Next, lets look at a change in the error
	foreach($newABCDError in $otherAbcdErrors) {
	    Start-LumiCalculation "$newABCDError ABCD Error" $inputFile $abcd $lumi $newABCDError $jobID $dataset
		#Write-Output "$dataset, $newABCDError ABCD Error, $($limit[0]), $($limit[1])"
	}

	foreach($abcdFactor in $abcdFactors) {
		$newabcd = Get-ABCDScaled $abcd $abcdFactor
	    Start-LumiCalculation "$abcdFactor background scaling" $inputFile $newabcd $lumi $abcdError $jobID $dataset
		#Write-Output "$dataset, $abcdFactor background scaling, $($limit[0]), $($limit[1])"
	}
}

# Finally, collect all the jobs and dump out the info.
$job | Wait-Job | Receive-Job | % { Write-Output "$_[1] $_[2] $_[0]" }
