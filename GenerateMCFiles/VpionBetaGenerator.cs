using libDataAccess;
using libDataAccess.Utils;
using LINQToTreeHelpers.FutureUtils;
using LinqToTTreeInterfacesLib;
using LINQToTTreeLib;
using LINQToTTreeLib.CodeAttributes;
using LINQToTTreeLib.Files;
using System;
using System.IO;
using System.Linq;
using static System.Math;
using System.Collections.Generic;

namespace GenerateMCFiles
{
    static class VpionBetaGenerator
    {
        /// <summary>
        /// Generate a ROOT file that will have all the MC events in them along with the vpion stuff.
        /// </summary>
        /// <param name="sample">The sample of data to skim</param>
        /// <param name="namePostfix">The output filename post-fix we should add.</param>
        public static IFutureValue<FileInfo> CreateVpionROOTFiles(this SampleMetaData sample, string namePostfix)
        {
            // Get the sample file
            var file = Files.GetSampleAsMetaData(sample);

            // The object that gives us access to the official event selection and ABCD method region
            // boundaries.
            var eventSelector = new SelectionHelperEventSelection();

            // Create the data we are going to write out. Every single event
            // has to be written out.
            // TODO:
            //   - Change to BDT13Lxy
            //   - Change to using BDT ordered jets
            //  Both of these require upgrading to the next version of the TTree, sadly.
            var dataStream = from evt in file
                             where evt.Data.LLPs.Count() == 2
                             let llp1 = evt.Data.LLPs[0]
                             let llp2 = evt.Data.LLPs[1]
                             let jets = evt.Data.Jets
                                            .OrderByDescending(j => j.BDT13Lxy)
                                            .Take(2)
                             let j1 = jets.FirstOrDefault()
                             let j2 = jets.Skip(1).FirstOrDefault()
                             let minDR2Sum = evt.Data.Jets
                                                .Where(j => j.pT > 50.0 && Abs(j.eta) < 2.5 && j.isGoodLLP)
                                                .Sum(j => j.CalibJet_minDRTrkpt2)
                             let isSelected = jets.Count() != 2
                                ? false
                                : eventSelector.eventSelection(j1.pT, j2.pT, j1.eta, j2.eta, j1.isGoodLLP, j2.isGoodLLP,
                                    j1.phi, j2.phi, j1.time, j2.time, evt.Data.event_HTMiss, evt.Data.event_HT, minDR2Sum)
                             let passedTrigger = evt.Data.event_passCalRatio_TAU60
                             let region = (jets.Count() != 2)
                                ? 0
                                : eventSelector.ABCDPlane(j1.pT, j2.pT, j1.BDT13Lxy, j2.BDT13Lxy, minDR2Sum)
                             select new VpionData
                             {
                                 PassedCalRatio = passedTrigger,
                                 eventNumber = evt.Data.eventNumber,
                                 llp1_E = llp1.E,
                                 llp1_eta = llp1.eta,
                                 llp1_phi = llp1.phi,
                                 llp1_pt = llp1.pT,
                                 llp1_Lxy = llp1.Lxy,
                                 llp2_E = llp2.E,
                                 llp2_eta = llp2.eta,
                                 llp2_phi = llp2.phi,
                                 llp2_pt = llp2.pT,
                                 llp2_Lxy = llp2.Lxy,
                                 event_weight = evt.Data.eventWeight * Abs(evt.Data.pileupEventWeight),
                                 // TODO: get from Emma how to do this correctly (once we figure it out!!)
                                 RegionA = passedTrigger && isSelected && region == 1,
                                 RegionB = passedTrigger && isSelected && region == 2,
                                 RegionC = passedTrigger && isSelected && region == 3,
                                 RegionD = passedTrigger && isSelected && region == 4
                             };

            // Now, write it out to a file.
            var f = dataStream
                .FutureAsTTree(treeName: "extrapTree", treeTitle: "Used as input for the extrapolation");

            // Total up everything for final numbers, and dump it out.
            var total = dataStream.FutureCount();
            var totalA = dataStream.Where(t => t.RegionA).FutureCount();
            var totalB = dataStream.Where(t => t.RegionB).FutureCount();
            var totalC = dataStream.Where(t => t.RegionC).FutureCount();
            var totalD = dataStream.Where(t => t.RegionD).FutureCount();

            FutureConsole.FutureWriteLine(() => $"EventInfo: {namePostfix} {total.Value} {totalA.Value} {totalB.Value} {totalC.Value} {totalD.Value}");

            // Total up everything for final numbers, and dump it out. Use weights
            var totalW = dataStream.FutureAggregate(0.0, (ac, ev) => ac + ev.event_weight);
            var totalAW = dataStream.Where(t => t.RegionA).FutureAggregate(0.0, (ac, ev) => ac + ev.event_weight);
            var totalBW = dataStream.Where(t => t.RegionB).FutureAggregate(0.0, (ac, ev) => ac + ev.event_weight);
            var totalCW = dataStream.Where(t => t.RegionC).FutureAggregate(0.0, (ac, ev) => ac + ev.event_weight);
            var totalDW = dataStream.Where(t => t.RegionD).FutureAggregate(0.0, (ac, ev) => ac + ev.event_weight);

            FutureConsole.FutureWriteLine(() => $"EventInfo (weighted): {namePostfix} {totalW.Value:F1} {totalAW.Value:F1} {totalBW.Value:F1} {totalCW.Value:F1} {totalDW.Value:F1}");

#if DEBUG_TEST
            var csvdata = dataStream
                .Where(e => e.RegionA || e.RegionB || e.RegionC || e.RegionD)
                .FutureAsCSV(new FileInfo("events-in-signal-region.csv"));
#endif

            // Return only the first file - as there should be no more than that!
            return f
                .Select(flst =>
                    {
                        if (flst.Length != 1)
                        {
                            throw new ArgumentException($"Got more than one file back when running LLP Extrapolation ntuple (found {flst.Length}!");
                        }
                        return flst;
                    })
                 .Select(flst => flst[0])
                 .Select(fspec =>
                 {
                     var fnew = new FileInfo(Path.Combine(fspec.DirectoryName, $"LLPExtrapolationMCTree-{namePostfix}{fspec.Extension}"));
                     return fspec.CopyTo(fnew.FullName, true);
                 });
        }

        /// <summary>
        /// Define the shape of the output ntuple.
        /// </summary>
        public class VpionData
        {
            public int eventNumber;
            /// <summary>
            /// Did it pass the CalRatio trigger?
            /// </summary>
            public bool PassedCalRatio;

            /// <summary>
            /// Info about the first LLP
            /// </summary>
            public double llp1_pt;
            public double llp1_eta;
            public double llp1_phi;
            public double llp1_E;
            public double llp1_Lxy;

            /// <summary>
            /// Info about the second LLP
            /// </summary>
            public double llp2_pt;
            public double llp2_eta;
            public double llp2_phi;
            public double llp2_E;
            public double llp2_Lxy;

            /// <summary>
            /// Any MC event wight
            /// </summary>
            /// <remarks>Seems to be 1.0 for this round of the analysis</remarks>
            public double event_weight;

            /// <summary>
            /// What ABCD region is it in? It is possible to be in no region at all (not even in the plane).
            /// </summary>
            public bool RegionA;
            public bool RegionB;
            public bool RegionC;
            public bool RegionD;
        }
    }

    /// <summary>
    /// Generate CPP code for the selection function and the ABCD region finder.
    /// </summary>
    /// <remarks>
    /// A regular CPP code static object would be just as good, but if
    /// we do this we can include the code directly as a comment, and then when it changes, force a re-run
    /// since the hash of the function text will change.
    /// </remarks>
    class SelectionHelperEventSelection : IOnTheFlyCPPObject
    {
        /// <summary>
        /// Location of the include file.
        /// </summary>
        readonly string _include_file = @"..\..\..\..\..\CalRSkimmer\CalRSelection.h";

        /// <summary>
        /// Return that file as an include file.
        /// </summary>
        /// <returns></returns>
        public string[] IncludeFiles()
        {
            return new[] { _include_file };
        }

        /// <summary>
        /// Return the lines of code. Include the file as a series of comments so they are hashed as part of the query.
        /// </summary>
        /// <param name="methodName"></param>
        /// <returns></returns>
        public IEnumerable<string> LinesOfCode(string methodName)
        {
            // First, generate the commented out code.
            using (var reader = File.OpenText(_include_file))
            {
                string line;
                while ((line = reader.ReadLine()) != null)
                {
                    yield return "// " + line;
                }
            }

            // Emit the specific code depending on which method is getting called.
            if (methodName == "eventSelection")
            {
                yield return "eventSelection = event_selection(j1_pt, j2_pt, j1_eta, j2_eta, j1_isGoodLLP, j2_isGoodLLP, j1_phi, j2_phi, j1_time, j2_time, event_HTMiss, event_HT,sumMinDRTrk2pt50);";
            }
            else if (methodName == "ABCDPlane")
            {
                yield return "ABCDPlane = event_ABCD_plane(j1_pt, j2_pt, j1_bdt13lxy, j2_bdt13lxy, sumMinDRTrk2pt50);";
            }
            else
            {
                throw new InvalidOperationException($"Unknown method {methodName} referenced - don't know how to code it!");
            }
        }

        /// <summary>
        /// The official base event selection.
        /// </summary>
        /// <param name="j1_pt"></param>
        /// <param name="j2_pt"></param>
        /// <param name="j1_eta"></param>
        /// <param name="j2_eta"></param>
        /// <param name="j1_isGoodLLP"></param>
        /// <param name="j2_isGoodLLP"></param>
        /// <param name="j1_phi"></param>
        /// <param name="j2_phi"></param>
        /// <param name="j1_time"></param>
        /// <param name="j2_time"></param>
        /// <param name="event_HTMiss"></param>
        /// <param name="event_HT"></param>
        /// <returns></returns>
        public bool eventSelection(
            double j1_pt, double j2_pt,
            double j1_eta, double j2_eta,
            bool j1_isGoodLLP, bool j2_isGoodLLP,
            double j1_phi, double j2_phi,
            double j1_time, double j2_time,
            double event_HTMiss, double event_HT,
            double sumMinDRTrk2pt50)
        {
            throw new InvalidOperationException("Should never get called by C# code!");
        }

        /// <summary>
        /// Which region of the ABCD plane is the event in?
        /// </summary>
        /// <param name="j1_pt"></param>
        /// <param name="j2_pt"></param>
        /// <param name="j1_bdt13lxy"></param>
        /// <param name="j2_bdt13lxy"></param>
        /// <param name="sumMinDRTrk2pt50"></param>
        /// <returns></returns>
        public int ABCDPlane(
            double j1_pt, double j2_pt,
            double j1_bdt13lxy, double j2_bdt13lxy,
            double sumMinDRTrk2pt50)
        {
            throw new InvalidOperationException("Should never get called by C# code!");
        }
    }
}
