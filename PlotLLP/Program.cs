using LINQToTreeHelpers.FutureUtils;
using LINQToTTreeLib;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static LINQToTreeHelpers.PlottingUtils;

namespace PlotLLP
{
    /// <summary>
    /// Generate a set of plots for the LLP's
    /// </summary>
    class Program
    {
        static void Main(string[] args)
        {
            // Get the input file from the command line.
            if (args.Length != 1)
            {
                throw new ArgumentException("Usage: PlotLLP <input-filename>");
            }
            var f = new FileInfo(args[0]);
            if (!f.Exists)
            {
                throw new ArgumentException($"Unable to find file {f.FullName}.");
            }

            // Get the query guy we want to plot
            var q = LLPInfo.QueryableextrapTree.CreateQueriable(f);

            // And the output file.
            var rootFilename = Path.GetFileNameWithoutExtension(f.Name);
            var outputFile = new FileInfo($"info-{rootFilename}.root");
            using (var of = new FutureTFile(outputFile))
            {
                // Do a few generic plots on the whole thing
                q
                    .GenericPlots(of, "all");

                // Right around the center where there are some real differences in Lxy
                q
                    .GenericPlotsInSquare(2000, 3500, of);
                q
                    .GenericPlotsInSquare(2500, 3000, of);
                q
                    .Where(evt => evt.RegionA)
                    .GenericPlots(of, "regionA");
                q
                    .Where(evt => evt.RegionA)
                    .GenericPlotsInSquare(2000, 3500, of, "A");
                q
                    .Where(evt => evt.RegionA)
                    .GenericPlotsInSquare(2500, 3000, of, "A");

            }
        }
    }

    static class PlotUtils
    {
        public static IPlotSpec<double> JetPtPlotRaw =
            MakePlotterSpec<double>(50, 0.0, 750.0, j => j, "pT{0}", "pT of {0} jets; pT [GeV]");

        public static IPlotSpec<double> JetLxyPlotRaw =
            MakePlotterSpec<double>(50, 0.0, 10, j => j, "Lxy{0}", "Lxy of {0} jets; Lxy [m]");

        public static IPlotSpec<double> JetEtaPlotRaw =
            MakePlotterSpec<double>(50, -6.0, 6.0, j => j, "Eta{0}", "Eta of {0} jets; Eta");

        public static void GenericPlotsInSquare(this IQueryable<LLPInfo.extrapTree> source, double lxy_low, double lxy_high, FutureTDirectory outd, string namearg = "")
        {
            source
                .Where(e => e.llp1_Lxy > lxy_low && e.llp1_Lxy< lxy_high)
                .Where(e => e.llp2_Lxy > lxy_low && e.llp2_Lxy< lxy_high)
                .GenericPlots(outd, $"{namearg}lxy_{lxy_low}_{lxy_high}");
    }

    /// <summary>
    /// Generate the generic plots
    /// </summary>
    /// <param name="source"></param>
    /// <param name="outd"></param>
    public static void GenericPlots(this IQueryable<LLPInfo.extrapTree> source, FutureTDirectory outd, string name)
        {
            source
                .Select(e => e.llp1_pt / 1000.0)
                .FuturePlot(JetPtPlotRaw, $"{name}_pt1", $"{name}_pt1")
                .Save(outd);
            source
                .Select(e => e.llp2_pt / 1000.0)
                .FuturePlot(JetPtPlotRaw, $"{name}_pt2", $"{name}_pt2")
                .Save(outd);
            source
                .Select(e => e.llp1_Lxy / 1000.0)
                .FuturePlot(JetLxyPlotRaw, $"{name}_lxy1", $"{name}_lxy1")
                .Save(outd);
            source
                .Select(e => e.llp2_Lxy / 1000.0)
                .FuturePlot(JetLxyPlotRaw, $"{name}_lxy2", $"{name}_lxy2")
                .Save(outd);
            source
                .Select(e => e.llp1_eta)
                .FuturePlot(JetEtaPlotRaw, $"{name}_eta1", $"{name}_eta1")
                .Save(outd);
            source
                .Select(e => e.llp2_eta)
                .FuturePlot(JetEtaPlotRaw, $"{name}_eta2", $"{name}_eta2")
                .Save(outd);
        }
    }
}
