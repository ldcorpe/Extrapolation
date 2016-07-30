using libDataAccess;
using LINQToTreeHelpers.FutureUtils;
using LinqToTTreeInterfacesLib;
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
        public static IFutureValue<FileInfo> CreateVpionROOTFiles (this SampleMetaData sample, string namePostfix)
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
                             let isSignal = jets.Count() != 2 
                                ? false
                                : j1.logRatio > 1.2 && j2.logRatio > 1.2 && j1.nTrk == 0 && j2.nTrk == 0
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
                                 RegionA = isSignal,
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
}
