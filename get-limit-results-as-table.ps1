#
# Grab the limits from jenksins along with parameters and emit an xml file.
#

# What range of jobs should we look at? Bad jobs will be skipped.

$jobMin = 18
$jobMax = 62

# Config - hopefully no need to modify

$jobUrl = "http://jenks-higgs.phys.washington.edu:8080/view/LLP/job/Limit-RunLimitExtrapolation/"

# Get the jobs and thus the parameters

$jobs = Find-JenkinsJob -JobUri $jobUrl -MinimumJobNumber $jobMin -MaximumJobNumber $jobMax -JobState Success

Write-Output "Job ID `t Rescale Signal `t Number of Toys `t 95% CL Rescaled `t 95% CL Raw `t Asym Fit"
foreach ($j in $jobs) {

	# Get the scaled results for each job
	$logFile = Get-JenkinsBuildLogfile $j
	$rescaled = $logFile | ? {$_.Contains("Limit rescaled: Obs:")}
	$raw = $logFile | where {$_.Contains("-> 95%=")}

	$cl95RS = -1.0
	if ($rescaled -match "95%: ([0-9\.]+)") {
		$cl95RS = $Matches[1]
	}

	$cl59Raw = -1.0
	if ($raw -match "-> 95%=([0-9\.]+)") {
		$cl59Raw = $Matches[1]
	}

	# Everything we want to write out
	$outJobId = $j.Id
	$outEventScaling = $j.Parameters["RescaleSignal"]
	$outIsAsym = $j.Parameters["UseAsymFit"]
	$outToys = 2500
	if ($j.Parameters.ContainsKey("NToys")) {
		$outToys = $j.Parameters["NToys"]
	}
	
	Write-Output "$outJobId `t $outEventScaling `t $outToys `t $cl95RS `t $cl59Raw `t $outIsAsym"	
}