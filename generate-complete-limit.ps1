#
# Given an MC files job, generate all the files downstream that are needed for limit plots.
#

# Base everything on this set of MC files
$mcFilesJobID = 32

# When we run a limit extrapolation, we need a template job. This is the id of the template job which we will
# base all our runs off of.
$runLimitExtrapTemplateJobID = 73

# Minimum job number in the LifetimeEfficiencies and limit extrap jenkins task to search, just to keep things reasonable
$mcEfficienciesMinJobID = 95
$runLimitExtrapMinJobID = 70

#
# Other config (hopefully don't need to be modified much!)
#

$mcFilesJobUri = "http://jenks-higgs.phys.washington.edu:8080/view/LLP/job/Limit-MCFiles/"
$mcEfficienciesJobUri = "http://jenks-higgs.phys.washington.edu:8080/view/LLP/job/Limit-LifetimeEfficiencies/"
$mcRunLimitExtrap = "http://jenks-higgs.phys.washington.edu:8080/view/LLP/job/Limit-RunLimitExtrapolation/"
$mcInfo = Import-Csv "$PSScriptRoot\GenerateMCFiles\Sample Meta Data.csv"

#
# Processing
#

# Get the list of artifacts that came out of the MCFiles. These will have the limit guys that will
# drive everything from here on out.

$j = Find-JenkinsJob -JobUri $mcFilesJobUri -JobId $mcFilesJobID
if ($j.Status -ne "SUCCESS") {
	throw "Job $mcFilesJobID did not finish successfully!"
}
$artifacts = $j.Artifacts | % { $_.ToString().TrimEnd(".root").TrimStart("LLPExtrapolationMCTree-") }

# Helper function that will look for a job for a dataset that has the efficiencies done.
# It will submit them if it hasn't been done.
$jobsMCEff = Find-JenkinsJob -JobUri $mcEfficienciesJobUri -MinimumJobNumber $mcEfficienciesMinJobID -ExcludeRebuiltJobs
function findMCEff ($dataset)
{
	# Find a matching job
	$jobMatch = $jobsMCEff | where {$_.Parameters["JobID"] -eq $mcFilesJobID} | where {$_.Parameters["Dataset"] -eq $dataset}

	# If we have no job, then submit it.
	if (!$jobMatch) {
		$p = @{}
		$p["Dataset"] = $dataset.ToString()
		$p["JobID"] = "$mcFilesJobID"
		$info = $mcInfo | where {$_."Nick Name" -eq $dataset}
		if (!$info) {
			throw "MC Dataset with nickname $a is not in the master csv file - can't figure out its lifetime!"
		}
		$pl = $info."Proper Lifetime [m]".ToString()
		$p["ProperLifetime"] = $pl
		Invoke-JenkinsJob -JobUri $mcEfficienciesJobUri -ParameterValues $p
		return "$dataset eff projection was submitted"
	} else {
		if ($jobMatch.Status -eq "SUCCESS") {
			return $jobMatch.Id
		} else {
			return "MC Efficiencies job $($jobMatch.Id) for $dataset is not completed yet ($($jobMatch.Status)). Either wait or rebuild with fix"
		}
	}
}

# Helper function that will look for a job in the Limit calculation, starting from the job id
# of the efficiency extrapolation
$jobsLimitExtrap = Find-JenkinsJob -JobUri $mcRunLimitExtrap -MinimumJobNumber $runLimitExtrapMinJobID -ExcludeRebuiltJobs
function findLimitExtrap($jobid)
{
	# If we are dealing with something that isn't an int, just pass it through.
	if ($jobid -notmatch "^[\d\.]+$") {
		return $jobid
	}

	# Now, see if we can find a job.
	$jobMatch = $jobsLimitExtrap | where {$_.Parameters["JobID"] -eq $jobid}

	$origJob = Find-JenkinsJob -JobUri $mcEfficienciesJobUri -JobId $jobid
	$ds = $origJob.Parameters["Dataset"]
	if (!$jobMatch) {
		$templateJob = Find-JenkinsJob -JobUri $mcRunLimitExtrap -JobId $runLimitExtrapTemplateJobID
		$p = $templateJob.Parameters
		$p["JobID"] = $jobid.ToString()
		$p["Dataset"] = $ds
		Invoke-JenkinsJob -JobUri $mcRunLimitExtrap -ParameterValues $p
		return "Limit Extrapolation job for $ds submitted."
	}
	if ($jobMatch.Status -ne "SUCCESS") {
		return "Limit Extrapolation job for $ds ($($jobMatch.Id)) has not finished ($($jobMatch.Status)). Please re-run script or fix job error."
	}
	$r = @{}
	$r["Dataset"] = $ds
	$r["Id"] = $jobMatch.Id
	$r["EffJobId"] = $jobid
	return $r
}

# Pull the event counts from a line like "EventInfo (weighted): 200pi50lt9mW 188988.0 34.0 6.0 5.0 0.0"
function Get-EventCount ($line)
{
	$m = $line -match " (\w+) ([\.0-9]+) ([\.0-9]+) ([\.0-9]+) ([\.0-9]+) ([\.0-9]+)$"
	$p = @{"ds" = $Matches[1]; "total" = $Matches[2]; "A" = $Matches[3] };
	return $p
}

# Given the results from the final file, generate everything plot lingo needs to drive it.
function generatePlotLingoSnippit ()
{
	Process {
		if ($_ -is [string]) {
			return $_
		}

		$ds = $_.Dataset
		$jobid = $_.Id
		$effJobId = $_.EffJobId

		# First, the line that will open the file.
		Write-Output "f_limit_$ds = jenkins(""http://jenks-higgs.phys.washington.edu:8080/view/LLP/job/Limit-RunLimitExtrapolation/$jobid/artifact/limit_$ds.root"");"
		Write-Output "f_eff_$ds = jenkins(""http://jenks-higgs.phys.washington.edu:8080/view/LLP/job/Limit-LifetimeEfficiencies/$effJobId/artifact/lifetime_eff_$ds.root"");"

		# Add it to the list of datasets
		Write-Output "ds_list = ds_list + [""$ds""];"
	}
	End {
	
		# Grab the raw numbers of events from the input jobs
		$log = Get-JenkinsBuildLogfile -JobUri $mcFilesJobUri -JobId $mcFilesJobID
		$unweighted = $log | where {$_ -match "EventInfo:"}
		$weighted = $log | where {$_ -match "EventInfo \("}
		$unweighted | % {Get-EventCount $_} | % { Write-Output "evt_$($_.ds)_tot = $($_.total);"; Write-Output "evt_$($_.ds)_A = $($_.A);"; $ds = $_.ds; $info = $mcInfo | where {$_."Nick Name" -eq $ds}; $pl = $info."Proper Lifetime [m]".ToString(); Write-Output "ctau_$($_.ds) = $pl;"}
		$weighted | % {Get-EventCount $_} | % { Write-Output "evtw_$($_.ds)_tot = $($_.total);"; Write-Output "evtw_$($_.ds)_A = $($_.A);"; $ds = $_.ds; $info = $mcInfo | where {$_."Nick Name" -eq $ds}; $pl = $info."Proper Lifetime [m]".ToString(); Write-Output "ctauw_$($_.ds) = $pl;" }
	}
}

# Run it all as a pipe-line.
$joblist = $artifacts `
	| % {findMCEff($_)}  `
	| % {findLimitExtrap($_)} `
	| generatePlotLingoSnippit
Write-Output $joblist
