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
                                 vpi1_E = llp1.LLP_E,
                                 vpi1_eta = llp1.eta,
                                 vpi1_phi = llp1.phi,
                                 vpi1_pt = llp1.pT,
                                 vpi1_Lxy = llp1.Lxy,
                                 vpi2_E = llp2.LLP_E,
                                 vpi2_eta = llp2.eta,
                                 vpi2_phi = llp2.phi,
                                 vpi2_pt = llp2.pT,
                                 vpi2_Lxy = llp2.Lxy,
                                 event_weight = evt.Data.eventWeight,
                                 // TODO: get from Emma how to do this correctly (once we figure it out!!)
                                 IsInSignalRegion = isSignal
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
            public bool PassedCalRatio;
            public double vpi1_pt;
            public double vpi1_eta;
            public double vpi1_phi;
            public double vpi1_E;
            public double vpi1_Lxy;
            public double vpi2_pt;
            public double vpi2_eta;
            public double vpi2_phi;
            public double vpi2_E;
            public double vpi2_Lxy;
            public double event_weight;
            public bool IsInSignalRegion;
        }

    }
}
