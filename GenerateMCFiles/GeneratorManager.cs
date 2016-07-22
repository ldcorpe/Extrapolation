using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using static System.Console;

namespace GenerateMCFiles
{
    class GeneratorManager
    {
        /// <summary>
        /// Top level generator for a dataset
        /// </summary>
        /// <param name="ds"></param>
        public static void GenerateExtrapolationMCFiles(string ds)
        {
            // Make sure this is a real dataset name that we can identify.
            var dsSpec = libDataAccess.SampleMetaData.LoadFromCSV(ds);
            if (dsSpec == null)
            {
                throw new ArgumentException($"Dataset '{ds}' was not found in the official DS list!");
            }

            // Generate the ROOT file(s) that contain the basic vpion information for every event that is needed
            // to create the extrapolation.

            var f = dsSpec.CreateVpionROOTFiles(ds);
            WriteLine(f.Value.FullName);
        }
    }
}
