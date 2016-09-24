#
# This will run the official limit several times with different configurations.
#

#
# Config the default settings
#

# The number of events in each of the A, B, C, and D regions
$abcd = (4, 16, 39, 34)

# Lumi in fb of the current results
$lumi = 3.2
$otherLumis = (30, 300)

# the size of the ABCD error
$abcdError = 0.36
$otherAbcdErrors = (0.10, 0.20)

# Fraction of background we find.
$abcdFactors = (0.90, 0.80, 0.60)

# What samples to look at and where to get the MC root files.
$singalSamples = (("1000pi400lt5m", 104), ("400pi50lt5m", 116))

function Invoke-LumiRun ($abcdInfo, $lumi, $abcdError, $jobID, $dataset) {
	$r = .\Release\ExtrapLimitFinder.exe -extrapFile $inputFile.FullName -UseAsym -nA $abcd[0] -nB $abcd[1] -nC $abcd[2] -nD $abcd[3] -OutputFile $dataset-normal.root -Luminosity $lumi -ABCDError $abcdError
	#$r = Get-Content junk.txt
	$l = $r -match "Limit rescaled"
	if ($l.Trim() -match "95%: ([0-9]+\.[0-9]+)") {
		$limit = $Matches[1]
		return $limit
	}
	return 0.0
}

function Get-ABCDScaled ($abcdOriginal, $scale) {
	$a = $abcdOriginal[0] * $scale
	$b = $abcdOriginal[1] * $scale
	$c = $abcdOriginal[2] * $scale
	$d = $abcdOriginal[3] * $scale
	return ($a, $b, $c, $d)
}

#########
# Run it!

foreach ($signal in $singalSamples) {
	$dataset = $signal[0]
	$jobID = $signal[1]

	# Get the input file from jenkins
	$inputFile = Get-JenkinsArtifact -JobUri http://jenks-higgs.phys.washington.edu:8080/view/LLP/job/Limit-LifetimeEfficiencies/ -JobId $jobID -ArtifactName lifetime_eff_$dataset.root

	# Get the default running result
	$limit = Invoke-LumiRun $abcd $lumi $abcdError $jobID $dataset
	Write-Output "$dataset, Default, $limit"

	# Next, do the increased amount of lumi
	foreach ($newLumi in $otherLumis) {
		$factor = $newLimi / $lumi
		$newabcd = Get-ABCDScaled $abcd $factor
	    $limit = Invoke-LumiRun $newabcd $newLimi $abcdError $jobID $dataset
		Write-Output "$dataset, $newLumi fb, $limit"
	}

	# Next, lets look at a change in the error
	foreach($newABCDError in $otherAbcdErrors) {
	    $limit = Invoke-LumiRun $abcd $lumi $newABCDError $jobID $dataset
		Write-Output "$dataset, $newABCDError ABCD Error, $limit"
	}

	foreach($abcdFactor in $abcdFactors) {
		$newabcd = Get-ABCDScaled $abcd $abcdFactor
	    $limit = Invoke-LumiRun $newabcd $newLimi $abcdError $jobID $dataset
		Write-Output "$dataset, $abcdFactor background scaling, $limit"
	}
}
