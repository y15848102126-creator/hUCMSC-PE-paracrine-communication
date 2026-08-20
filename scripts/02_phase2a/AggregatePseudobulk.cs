using System;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Text;

public static class AggregatePseudobulk
{
    public static int Main(string[] args)
    {
        if (args.Length != 5)
        {
            Console.Error.WriteLine("Usage: AggregatePseudobulk ZIP GROUP_IDS STRATA_TSV OUTPUT_GZ TOTALS_TSV");
            return 2;
        }
        int[] groupIds;
        using (var reader = new BinaryReader(File.OpenRead(args[1])))
        {
            int n = checked((int)(reader.BaseStream.Length / 4));
            groupIds = new int[n];
            for (int i = 0; i < n; i++) groupIds[i] = reader.ReadInt32();
        }
        string[] strataLines = File.ReadAllLines(args[2], Encoding.UTF8);
        int groupN = strataLines.Length - 1;
        string[] columns = new string[groupN];
        for (int i = 0; i < groupN; i++) columns[i] = strataLines[i + 1].Split('\t')[3];
        long[] grandTotals = new long[groupN];
        long geneN = 0;
        var started = DateTime.UtcNow;
        using (var zip = ZipFile.OpenRead(args[0]))
        using (var input = new StreamReader(zip.Entries[0].Open(), Encoding.UTF8, true, 1 << 20))
        using (var outputFile = File.Create(args[3]))
        using (var gzip = new GZipStream(outputFile, CompressionLevel.Optimal))
        using (var output = new StreamWriter(gzip, new UTF8Encoding(false), 1 << 20))
        {
            for (int i = 0; i < 23; i++) input.ReadLine();
            output.Write("gene");
            for (int i = 0; i < groupN; i++) { output.Write(','); output.Write(columns[i]); }
            output.WriteLine();
            string line;
            while ((line = input.ReadLine()) != null)
            {
                int firstTab = line.IndexOf('\t');
                if (firstTab < 1) throw new InvalidDataException("Expression row lacks a tab");
                long[] sums = new long[groupN];
                int cell = 0;
                long value = 0;
                for (int pos = firstTab + 1; pos <= line.Length; pos++)
                {
                    char ch = pos == line.Length ? '\t' : line[pos];
                    if (ch == '\t')
                    {
                        if (cell >= groupIds.Length) throw new InvalidDataException("Too many cell values");
                        sums[groupIds[cell]] += value;
                        cell++;
                        value = 0;
                    }
                    else if (ch >= '0' && ch <= '9') value = checked(value * 10 + (ch - '0'));
                    else throw new InvalidDataException("Non-integer UMI value at gene " + line.Substring(0, firstTab));
                }
                if (cell != groupIds.Length) throw new InvalidDataException("Wrong cell count for gene " + line.Substring(0, firstTab));
                output.Write(line.Substring(0, firstTab).Replace("\"", "\"\""));
                for (int i = 0; i < groupN; i++)
                {
                    output.Write(','); output.Write(sums[i].ToString(CultureInfo.InvariantCulture)); grandTotals[i] += sums[i];
                }
                output.WriteLine();
                geneN++;
                if (geneN % 5000 == 0) Console.WriteLine("aggregated " + geneN + " genes in " + (DateTime.UtcNow - started).TotalSeconds.ToString("F1", CultureInfo.InvariantCulture) + "s");
            }
        }
        using (var totals = new StreamWriter(args[4], false, new UTF8Encoding(false)))
        {
            totals.WriteLine("group_id\tmatrix_column\tobserved_total_umi");
            for (int i = 0; i < groupN; i++) totals.WriteLine(i + "\t" + columns[i] + "\t" + grandTotals[i]);
        }
        Console.WriteLine("completed " + geneN + " genes in " + (DateTime.UtcNow - started).TotalSeconds.ToString("F1", CultureInfo.InvariantCulture) + "s");
        return 0;
    }
}
