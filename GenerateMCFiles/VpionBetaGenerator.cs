using libDataAccess;
using System;
using System.Collections.Generic;
using System.Linq;
using LINQToTTreeLib.Files;
using System.Text;
using System.Threading.Tasks;
using LinqToTTreeInterfacesLib;
using LINQToTreeHelpers.FutureUtils;
using LINQToTTreeLib;
using System.IO;

namespace GenerateMCFiles
{
    static class VpionBetaGenerator
    {
        /// <summary>
        /// Generate a ROOT file that will have all the MC events in them along with the vpion stuff.
        /// </summary>
        /// <param name="sample"></param>
        public static IFutureValue<FileInfo> CreateVpionROOTFiles (this SampleMetaData sample)
        {
            // Get the sample file
            var file = Files.GetSampleAsMetaData(sample);

            // Create the data we are going to write out. Every single event
            // has to be written out.
            var dataStream = from evt in file
                             where evt.Data.LLPs.Count() == 2
                             let llp1 = evt.Data.LLPs[0]
                             let llp2 = evt.Data.LLPs[1]
                             select new VpionData
                             {
                                 // TODO: Make sure from Rachel and Emma this is the right trigger - seems to be zero all the time.
                                 PassedCalRatio = evt.Data.event_passCalRatio_cleanLLP_TAU60,
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
                                 weight = 1.0,
                                 // TODO: get from Emma how to do this correctly (once we figure it out!!)
                                 IsInSignalRegion = true
                             };

            // Now, write it out to a file.
            var f = dataStream
                .FutureAsTTree(treeName: "extrapTree", treeTitle: "Used as input for the extrapolation");

            // Return only the first file - as there should be no more than that!
            return from flst in f select flst.First();
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
            public double weight;
            public bool IsInSignalRegion;
        }

    }
}
