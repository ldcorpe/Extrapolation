using libDataAccess;
using LINQToTreeHelpers.FutureUtils;
using LinqToTTreeInterfacesLib;
using LINQToTTreeLib.CodeAttributes;
using LINQToTTreeLib.Files;
using System;
using System.IO;
using System.Linq;
using static System.Math;

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

            // Create the data we are going to write out. Every single event
            // has to be written out.
            var dataStream = from evt in file
                             where evt.Data.LLPs.Count() == 2
                             let llp1 = evt.Data.LLPs[0]
                             let llp2 = evt.Data.LLPs[1]
                             let jets = evt.Data.Jets
                                            .Where(j => j.pT > 100.0 && Abs(j.eta) < 2.5)
                                            .OrderByDescending(j => j.pT)
                                            .Take(2)
                             let j1 = jets.FirstOrDefault()
                             let j2 = jets.Skip(1).FirstOrDefault()
                             let isSelected = jets.Count() != 2
                                ? false
                                : SelectionHelpers.eventSelection(j1.pT, j2.pT, j1.eta, j2.eta, j1.isGoodLLP, j2.isGoodLLP,
                                    j1.phi, j2.phi, j1.CalibJet_time, j2.CalibJet_time, evt.Data.event_HTMiss, evt.Data.event_HT)
                             select new VpionData
                             {
                                 PassedCalRatio = evt.Data.event_passCalRatio_TAU60,
                                 llp1_E = llp1.LLP_E,
                                 llp1_eta = llp1.eta,
                                 llp1_phi = llp1.phi,
                                 llp1_pt = llp1.pT,
                                 llp1_Lxy = llp1.Lxy,
                                 llp2_E = llp2.LLP_E,
                                 llp2_eta = llp2.eta,
                                 llp2_phi = llp2.phi,
                                 llp2_pt = llp2.pT,
                                 vpi2_Lxy = llp2.Lxy,
                                 event_weight = evt.Data.eventWeight,
                                 // TODO: get from Emma how to do this correctly (once we figure it out!!)
                                 RegionA = isSelected,
                                 RegionB = false,
                                 RegionC = false,
                                 RegionD = false
                             };

            // Now, write it out to a file.
            var f = dataStream
                .FutureAsTTree(treeName: "extrapTree", treeTitle: "Used as input for the extrapolation");

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
            public double vpi2_Lxy;

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
    /// There is C++ code that we share with everyone else. Use it here to make sure that
    /// the official selection cuts are identical everywhere in the analyiss.
    /// </summary>
    [CPPHelperClass()]
    class SelectionHelpers
    {
        /// <summary>
        /// Determine if the event is in the analysis (or not). J1 and J2 ordering are by BDT ordering, not
        /// by pT ordering as is standard in most analyses.
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
        [CPPCode(IncludeFiles = new[] { @"C:\Users\gordo\Documents\Code\calratio2015\CalRSkimmer\CalRSelection.h" }, Code = new[]
        {
            "eventSelection = event_selection(j1_pt, j2_pt, j1_eta, j2_eta, j1_isGoodLLP, j2_isGoodLLP, j1_phi, j2_phi, j1_time, j2_time, event_HTMiss, event_HT);"
        })]
        public static bool eventSelection(
            double j1_pt, double j2_pt,
            double j1_eta, double j2_eta,
            bool j1_isGoodLLP, bool j2_isGoodLLP,
            double j1_phi, double j2_phi,
            double j1_time, double j2_time,
            double event_HTMiss, double event_HT)
        {
            throw new InvalidOperationException("Should never get called by C# code!");
        }

        [CPPCode(IncludeFiles = new[] { "CalRSelection.h" }, Code = new[]
        {
            "ABCDPlane = event_ABCD_plane(j1_pt, j2_pt, j1_bdt13lxy, j2_bdt13lxy, sumMinDRTrk2pt50);"
        })]
        public static int ABCDPlane(
            double j1_pt, double j2_pt,
            double j1_bdt13lxy, double j2_bdt13lxy,
            double sumMinDRTrk2pt50)
        {
            throw new InvalidOperationException("Should never get called by C# code!");
        }
    }
}
