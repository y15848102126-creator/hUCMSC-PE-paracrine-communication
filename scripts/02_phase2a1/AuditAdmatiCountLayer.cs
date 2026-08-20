using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Text;

public static class AuditAdmatiCountLayer
{
    public static int Main(string[] args)
    {
        if (args.Length != 2) { Console.Error.WriteLine("Usage: AuditAdmatiCountLayer ZIP OUTPUT_CSV"); return 2; }
        string[] cellIds = null;
        long[] published = null;
        long[] observed = null;
        int[] detected = null;
        using (var zip = ZipFile.OpenRead(args[0]))
        using (var input = new StreamReader(zip.Entries[0].Open(), Encoding.UTF8, true, 1 << 20))
        {
            for (int row = 0; row < 23; row++)
            {
                string line = input.ReadLine();
                if (line == null) throw new InvalidDataException("Missing metadata rows");
                string[] fields = line.Split('\t');
                if (row == 0) { cellIds = new string[fields.Length - 1]; Array.Copy(fields, 1, cellIds, 0, cellIds.Length); observed = new long[cellIds.Length]; detected = new int[cellIds.Length]; }
                if (fields[0] == "total_molecules")
                {
                    published = new long[fields.Length - 1];
                    for (int i = 1; i < fields.Length; i++) published[i - 1] = long.Parse(fields[i], CultureInfo.InvariantCulture);
                }
            }
            if (published == null) throw new InvalidDataException("total_molecules row not found");
            string geneLine;
            long geneN = 0;
            while ((geneLine = input.ReadLine()) != null)
            {
                int firstTab = geneLine.IndexOf('\t');
                if (firstTab < 1) throw new InvalidDataException("Expression row lacks tab");
                int cell = 0; long value = 0;
                for (int pos = firstTab + 1; pos <= geneLine.Length; pos++)
                {
                    char ch = pos == geneLine.Length ? '\t' : geneLine[pos];
                    if (ch == '\t') { observed[cell] += value; if (value > 0) detected[cell]++; cell++; value = 0; }
                    else if (ch >= '0' && ch <= '9') value = checked(value * 10 + (ch - '0'));
                    else throw new InvalidDataException("Non-integer count");
                }
                if (cell != observed.Length) throw new InvalidDataException("Wrong cell count");
                geneN++;
                if (geneN % 5000 == 0) Console.WriteLine("audited " + geneN + " genes");
            }
        }
        using (var output = new StreamWriter(args[1], false, new UTF8Encoding(false)))
        {
            output.WriteLine("cell_id,published_total_molecules,direct_matrix_sum,detected_gene_n,direct_minus_detected_gene_n,consistent_with_ceil_to_10000,exact_raw_total");
            for (int i = 0; i < observed.Length; i++)
            {
                bool ceiling10k = observed[i] >= 10000 && observed[i] <= 10000 + detected[i];
                output.WriteLine(cellIds[i] + "," + published[i] + "," + observed[i] + "," + detected[i] + "," + (observed[i] - detected[i]) + "," + (ceiling10k ? "YES" : "NO") + "," + (published[i] == observed[i] ? "YES" : "NO"));
            }
        }
        Console.WriteLine("completed cells=" + observed.Length);
        return 0;
    }
}
