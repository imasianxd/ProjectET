using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using Jia;
using System.Threading;
using System.Diagnostics;
using System.ComponentModel;
using System.Windows.Forms; //hack for multi core threading.
using System.Drawing;
using System.Drawing.Imaging;
using System.Windows.Forms.DataVisualization.Charting;
using Newtonsoft.Json;
using OxyPlot;
using OxyPlot.Series;
using OfficeOpenXml;
using OfficeOpenXml.Style;

namespace ProjectET
{
	class findEpitope
	{
		static Library l;
		private string root, oldDirectory, curTime, resultDir;
		rawData[] deserializedDatas = new rawData[3];
		DebugSettings DS = new DebugSettings();

		public findEpitope()
		{
			//init
			l = new Library();
			l.setLogVerbose(1);
			root = Directory.GetCurrentDirectory();
			curTime = DateTime.Now.TimeOfDay.ToString();
			oldDirectory = root + "/Old_Runs/" + curTime;
			resultDir = root + "/Results";


			if (!File.Exists(root + "/Debug_Settings.json"))
			{
				DS.analysisTime = 0;
				DS.cutoff = 0.1;
				DS.isError = false;
				DS.stepsCompleted = -1;
				DS.nonPubSequence = true;
				writeDebug();
			}
			else
				readDebug();

			Stopwatch sw = new Stopwatch();
			sw.Start();

//			callEntryPeptide ();
//			List<parsemodified> p = parseModifiedResults ();
//			Environment.Exit (0);
			//  parseResult();
			// parseTabToJson();
			//  readJsonIntoMemory();
			// readAllDataIntoMemory();
			//   makeFigure();
			//  getPopulationData();
			//  getEpitopeParameter();
			//  outputResult();
			//   rankResult();
			//   generateReport();
			//   makeBexTable();
			//    compileProteinseq();

            //DS.stepsCompleted = 0;
            if (DS.isError)
            {
                Console.WriteLine("Last analysis did not terminate properly...attempt to recover the data and continue (y/n)?");
				string t = Convert.ToString (Console.ReadLine ().ToLower ());
                if (t == "y" || t == "yes")
                {
                    Console.WriteLine("Depending on how much epitopes was predicted, you may encounter am out of memory error");
                    Console.WriteLine("Now we don't want to re-run the analysis, so if you run into this error");
                    Console.WriteLine("please open Debug_Settings.json and modify the value of \"stepsCompleted\" to 4");
					if (DS.stepsCompleted > 7 && DS.stepsCompleted < 10) 
					{
						readAllDataIntoMemory ();
					}
					else if (DS.stepsCompleted > 4)
                    {
                        readJsonIntoMemory();
                    }
                }
                else
                {
                    Console.WriteLine("You choose not to recover, press any key to start fresh, the residual data of previous analysis will be in ./old_runs");
                    Console.WriteLine("If you changed your mind and want to recover, exit and restart the program. do not press any key.");
                    Console.ReadKey();
                    DS.stepsCompleted = -1;
                }
            }

            DS.isError = true;
            switch (DS.stepsCompleted)
            {
                case -1:
                    if (!checkFileStructure())
                        goto case 0;
                    else
                    {
                        Console.WriteLine("See above for directory errors, please fix");
                        l.writeLog("Check File Structure Failed");
                        break;
                    }
			case 0:
				Console.WriteLine ("Checking files...");
                    //Directory.CreateDirectory(oldDirectory);
					if (DS.nonPubSequence) 
					{
						//GetFiles on DirectoryInfo returns a FileInfo object.
						var Files = new DirectoryInfo(root+"/Protein_Sequence").GetFiles();
						//FileInfo has a Name property that only contains the filename part.
						List<string> filenames = new List<string> ();
						List<string> accessions = new List<string> ();
						for (int i = 0; i < Files.Count(); i++) {
							filenames.Add (Files [i].Name.Substring(0,Files[i].Name.IndexOf(".")) + ">Self_Sequence" + (i + 1).ToString ());
							accessions.Add (Files [i].Name.Substring(0,Files[i].Name.IndexOf(".")));
						}
						l.write (filenames, path: root + "/Protein_Accession_Dictionary.txt", delete: true);
						l.write (accessions, path: root + "/accession.txt", delete: true);
					}
					else
					{
	                    getSequence();
					}
                    DS.stepsCompleted = 1;
                    writeDebug();
                    goto case 1;
                case 1:
                    Console.WriteLine("Prediting Epitopes...");
                    predictEpitope();
                    DS.stepsCompleted = 2;
                    writeDebug();
                    goto case 2;
                case 2:
                    Console.WriteLine("Making sense of results...");
                    parseResult();
                    DS.stepsCompleted = 3;
                    writeDebug();
                    goto case 3;
                case 3:
                    Console.WriteLine("Writing to json...");
                    parseTabToJson();
                    DS.stepsCompleted = 4;
                    writeDebug();
                    goto case 4;
                case 4:
                    Console.WriteLine("Loading from json...");
                    readJsonIntoMemory();
                    DS.stepsCompleted = 5;
                    writeDebug();
                    goto case 5;
                case 5:
                    Console.WriteLine("Making graphs and heatmaps...");
                    makeFigure();
                    DS.stepsCompleted = 6;
                    writeDebug();
                    goto case 6;
                case 6:
                    Console.WriteLine("Calculating the population coverage...");
                    getPopulationData();
                    DS.stepsCompleted = 7;
                    writeDebug();
                    goto case 7;
                case 7:
                    Console.WriteLine("Predicting solubility...");
                    getEpitopeParameter();
                    DS.stepsCompleted = 8;
                    writeDebug();
                    goto case 8;
				case 8:
					Console.WriteLine ("Summarizing results...");
                    outputResult();
                    DS.stepsCompleted = 9;
                    writeDebug();
                    goto case 9;
                case 9:
                    Console.WriteLine("Outputting results...");
                    rankResult();
                    DS.stepsCompleted = 10;
                    writeDebug();
                    goto case 10;
                case 10:
                    Console.WriteLine("Generating report & colored tables...");
                    //generateReport();
                    makeBexTable();
                    DS.stepsCompleted = 11;
                    writeDebug();
                    goto case 11;
                case 11:
                    Console.WriteLine("Analysis complete!");
                    break;
                default:
                    l.writeLog("unexpected stepsCompleted, check debugSettings.json");
                    throw new Exception();
            }
            
			sw.Stop();
			DS.analysisTime = sw.Elapsed.TotalSeconds;
			DS.isError = false;
			DS.stepsCompleted = -1;
			writeDebug();
		}

		/// <summary>
		/// check to see if the file intergrity is intact before continue.
		/// </summary>
		/// <returns></returns>
		public bool checkFileStructure()
		{
			bool error = false;
			//getsequence() file structure
			if (!File.Exists(root + "/Protein_Accession_Dictionary.txt"))
			{ Console.WriteLine("Check to make sure accession.txt exists and is filled"); error = true; }

			//predcition() file structure
			if (!Directory.Exists(root + "/netMHCII-2.2"))
			{ Console.WriteLine("netMHCII-2.2 folder missing, ensure you extracted all of the files"); error = true; }
			if (!Directory.Exists(root + "/netMHCIIpan-3.0"))
			{ Console.WriteLine("netMHCIIpan-3.0 folder missing, ensure you extracted all of the files"); error = true; }
			if (!Directory.Exists(root + "/iedb"))
			{ Console.WriteLine("iedb folder missing, ensure you extracted all of the files"); error = true; }

			if (!Directory.Exists(root + "/Old_Runs"))
				Directory.CreateDirectory(root + "/Old_Runs");
			//Directory.CreateDirectory(oldDirectory);

			if (Directory.Exists(resultDir))
				Directory.Move(resultDir, oldDirectory);
			else
				Directory.CreateDirectory(oldDirectory);
			Directory.CreateDirectory(resultDir);

			return error;
		}

		/// <summary>
		/// get protein sequence from genbank and uniprot
		/// </summary>
		public void getSequence()
		{

			List<string> accessionsIn = l.read (path: root + "/Protein_Accession_Dictionary.txt");
			List<string> accessionOut = new List<string> ();
			foreach (string s in accessionsIn) {
				try{
				accessionOut.Add (s.Substring (0, s.IndexOf (">")));
				}
				catch (Exception e) {
					throw new Exception ("Check the protein_accession_dictionary.txt");
				}
			}
			l.write (accessionOut, path: root + "/accession.txt", delete:true);

			Console.WriteLine("Pulling sequence from online database...");
			l.writeLog("[INFO] grabbing protein sequence from database");

			if (File.Exists(root + "/accession.txt")) //if the accesion list is available.
			{
				List<string> input = l.read(path: root + "/accession.txt");
				if (input.Count < 1) { throw new Exception("nothing in accession.txt"); l.writeLog("nothing in accession.txt"); }

				if (Directory.Exists(root + "/Protein_Sequence"))
					Directory.Move(root + "/Protein_Sequence", oldDirectory + "/Protein_Sequence");
				Directory.CreateDirectory(root + "/Protein_Sequence");

				foreach (string s in input) //go through each accession
				{
					string url = "";
					int y;
					l.writeLog("[INFO]Getting data for accession: " + s);
					if (s.Length == 6 && !int.TryParse(s, out y)) //uniprot
						url = string.Format("http://www.uniprot.org/uniprot/{0}.fasta", s);
					else //gi
						url = string.Format("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={0}&rettype=fasta&retmode=fasta", s);//seq_start=1&seq_stop=825")
					//some random accession might break this.....                                                                                                                           //outpuzt.Add (l.getDataFromAPI (url));

					l.write(l.getDataFromAPI(url), delete: false, path: string.Format("{1}/Protein_Sequence/{0}.fasta", s, root));
				}
			}
			else
			{
				l.writeLog("[ERR]Protein accession list is not found at ./accession.txt");
				throw new Exception("accession list not found, make sure it exists"); //if list not found.
			}

			Console.WriteLine("Pulling sequence from online database...DONE");
		}

		/// <summary>
		/// following code is used to call netMHCII and iedb tools to begin epitope prediction
		/// </summary>
		#region Prediction
		public void predictEpitope()
		{
			BackgroundWorker v2 = new BackgroundWorker();//thread worker for netmhcii2.2 and iedb
			v2 = new BackgroundWorker();
			v2.WorkerReportsProgress = true;
			v2.WorkerSupportsCancellation = true;
			Control.CheckForIllegalCrossThreadCalls = false;

			v2.DoWork += new DoWorkEventHandler(do2);

			if (v2.IsBusy != true)
				v2.RunWorkerAsync();
			else
			{
				l.writeLog("[ERR]Thread Busy");
				throw new Exception("Thread Busy...Rerun script. If it still give same error. restart your computer");
			}

			BackgroundWorker v3 = new BackgroundWorker(); //thread worker for netmhcii3.0
			v3 = new BackgroundWorker();
			v3.WorkerReportsProgress = true;
			v3.WorkerSupportsCancellation = true;
			Control.CheckForIllegalCrossThreadCalls = false;

			v3.DoWork += new DoWorkEventHandler(do3);

			if (v3.IsBusy != true)
				v3.RunWorkerAsync();
			else
			{
				l.writeLog("[ERR]Thread Busy");
				throw new Exception("Thread Busy...Rerun script. If it still give same error. restart your computer");
			}

			while (v2.IsBusy || v3.IsBusy) //hang the threads during prediction
				Thread.Sleep(10000);
		}
		private void do2(object sender, DoWorkEventArgs e)
		{
			callEntry(root + "/netMHCII-2.2");
			callEntry(root + "/iedb");
		}
		private void do3(object sender, DoWorkEventArgs e)
		{
			string dir = Directory.GetCurrentDirectory();
			callEntry(root + "/netMHCIIpan-3.0");
		}
		public void callEntry(string dir) //function used to begin the actual prediction
		{
			Console.WriteLine("Starting Prediction...");
			if (File.Exists(dir + "/entry"))
			{
				//timer to count the analysis time
				Stopwatch timer = new Stopwatch();
				timer.Start();

				if (Directory.Exists(dir + "/results")) //if data from an old run exists, move it
				{
					Directory.Move(dir + "/results", oldDirectory + "/" + dir.Substring(dir.LastIndexOf("/") + 1) + "_Results");
					Directory.CreateDirectory(dir + "/results");
					l.writeLog("Moved existing prediction result to " + oldDirectory + "/" + dir.Substring(dir.LastIndexOf("/") + 1) + "_Results");
				}

				//beging a new bash process and call the entry script
				Process p = new Process();
				p.StartInfo.WorkingDirectory = dir;
				p.StartInfo.Arguments = "entry";
				p.StartInfo.FileName = "bash";
				p.StartInfo.UseShellExecute = false;
				p.StartInfo.RedirectStandardOutput = true;
				p.StartInfo.RedirectStandardError = true;
				p.OutputDataReceived += (s, err) => Console.WriteLine(err.Data);
				p.ErrorDataReceived += (s, err) => Console.WriteLine(err.Data);
				p.Start();
				p.BeginOutputReadLine();
				p.BeginErrorReadLine();
				p.WaitForExit();

				timer.Stop();//stop timer and print time elapsed
				Console.WriteLine("Starting Prediction... DONE");
				Console.WriteLine("This analysis took a total of {0} Seconds", timer.Elapsed.ToString());
			}
			else
			{
				l.writeLog("[ERR]Unable to find the Entry Script in " + dir);
				throw new Exception("entry script does not exist, check log");
			}
		}

		public void callEntryPeptide() //function used to begin the actual prediction
		{
			string dir = root + "/netMHCII-2.2";
			Console.WriteLine("Starting Prediction...");
			if (File.Exists(dir + "/entry.peptide"))
			{
				//timer to count the analysis time
				Stopwatch timer = new Stopwatch();
				timer.Start();

				if (Directory.Exists(dir + "/modified_results")) //if data from an old run exists, move it
				{
					//Directory.Move(dir + "/modified_results/", oldDirectory+"/modified_results/");
					Directory.Delete (dir + "/modified_results", true);
					Directory.CreateDirectory(dir + "/modified_results");
					l.writeLog("deleted!! existing modified prediction result to " + oldDirectory + "/" + dir.Substring(dir.LastIndexOf("/") + 1) + "_modified_Results");
				}
				else
					Directory.CreateDirectory(dir + "/modified_results");

				//beging a new bash process and call the entry script
				Process p = new Process();
				p.StartInfo.WorkingDirectory = dir;
				p.StartInfo.Arguments = "entry.peptide";
				p.StartInfo.FileName = "bash";
				p.StartInfo.UseShellExecute = false;
				p.StartInfo.RedirectStandardOutput = true;
				p.StartInfo.RedirectStandardError = true;
				p.OutputDataReceived += (s, err) => Console.WriteLine(err.Data);
				p.ErrorDataReceived += (s, err) => Console.WriteLine(err.Data);
				p.Start();
				p.BeginOutputReadLine();
				p.BeginErrorReadLine();
				p.WaitForExit();

				timer.Stop();//stop timer and print time elapsed
				Console.WriteLine("Starting Prediction... DONE");
				Console.WriteLine("This analysis took a total of {0} Seconds", timer.Elapsed.ToString());
			}
			else
			{
				l.writeLog("[ERR]Unable to find the Entry Script in " + dir);
				throw new Exception("entry script does not exist, check log");
			}
		}
		#endregion

		/// <summary>
		/// following code is used to parse the individual prediction results into one big .tab file that can be opened in excel
		/// </summary>
		#region Covert Raw Results to TAB then To json
		public void parseResult() //after prediction
		{
			examine32BitWorkaround(root + "/netMHCII-2.2/results/", "netMHCII-2.2");
			examine32BitWorkaround(root + "/netMHCIIpan-3.0/results/", "netMHCIIpan-3.0");
			examine32BitWorkaround(root + "/iedb/results/", "iedb");
		}

		//we are going to rewrite this function...
		public void examine32BitWorkaround(string workingDir, string version) //build for low RAM use, so many GC.collect everywhere
		{
			Console.WriteLine("starting result summarization for version {0}...", version);
			List<string> output = new List<string>();


			Dictionary<string, string> accessionToAllergen = new Dictionary<string, string>();
			if (File.Exists(root + "/Protein_Accession_Dictionary.txt"))
			{
				foreach (string s in l.read(path: root + "/Protein_Accession_Dictionary.txt"))
				{
					if (!s.Contains(">") && !s.Contains("#") && s.Length < 2)
					{
						l.writeLog("[ERR]Improperly formated Protein_Accession_Dictionary.txt");
						throw new Exception("Improperly formated Protein_Accession_Dictionary.txt");
					}
					else if (s.Contains("#") || s.Length < 2)
					{
						//ignore comment and empty lines
					}
					else
						accessionToAllergen.Add(s.Substring(0, s.IndexOf(">")), s.Substring(s.IndexOf(">") + 1));
				}
				if (accessionToAllergen.Count == 0)
					l.writeLog("[WARN]nothing in Protein_Accession_Dictionary.txt");
			}
			else
				l.writeLog("[INFO]protein accession dictionary not found");

			string protein;

			foreach (string n in Directory.GetDirectories(workingDir)) //this step loads everything into memory
			{

				try
				{
					protein = (accessionToAllergen[n.Substring(n.LastIndexOf("/") + 1).Substring(0, n.Substring(n.LastIndexOf("/") + 1).IndexOf("."))]);
				}
				catch (Exception e)
				{
					protein = n.Substring(n.LastIndexOf("/") + 1);
					l.writeLog("[WARN]Accession was not found in dictionary, file name used for identifier. Accession: " + protein);
				}

				Dictionary<string, List<string[]>> alleles = new Dictionary<string, List<string[]>>();

				foreach (string f in Directory.GetFiles(n))
				{
					List<string> input = l.read(path: f);
					List<string> filter = new List<string>();
					bool start = false;
					string allele = f.Substring(f.LastIndexOf("/") + 1);
					Console.WriteLine("Processing: " + protein + "/" + allele);
					if (input.Count == 0)
						throw new Exception("error 0j01: " + protein + "/" + allele);

					if (version != "iedb")
					{
						for (int i = 0; i < input.Count; i++)
						{
							if (input[i].Contains("----------------") && start)
								break;

							if (input[i].Contains("----------------"))
							{
								i += 3;
								start = true;
							}

							if (start)
								filter.Add(input[i].Replace(" ", ","));
						}
					}
					else
						for (int i = 1; i < input.Count; i++)
					{
						filter.Add(input[i].Replace("\t", ","));
					}
					if (filter.Count == 0)
						throw new Exception("error 0j02: " + protein + "/" + allele);

					List<string[]> filter2 = new List<string[]>();

					foreach (string s in filter)
					{
						string[] temp = s.Split(char.Parse(","));
						temp = temp.Where(w => w != "").ToArray(); //data for 1 epitope
						filter2.Add(temp); //list of data for all epitope for 1 allele
					}
					alleles.Add(allele, filter2); //allele - all epitope
				}

				string tmp = "";

				output.Add("Antigen: " + protein); //print protein

				Dictionary<string, List<string>> epitopes_bindingcoreDic = new Dictionary<string, List<string>>();
				//foreach (KeyValuePair<string, List<string[]>> a in alleles) //print epitope
				for (int a = 0; a < alleles.Count; a++)
				{
					tmp = "\t"; //display allele
					//allele
					for (int i = 0; i < alleles.ElementAt(a).Value.Count; i++) //all epitope for 1 allele
					{

						if (version == "netMHCII-2.2")
						{
							if (!epitopes_bindingcoreDic.ContainsKey(alleles.ElementAt(a).Value[i][2]))
							{
								epitopes_bindingcoreDic.Add(alleles.ElementAt(a).Value[i][2], new List<string>());
								epitopes_bindingcoreDic[alleles.ElementAt(a).Value[i][2]].Add(alleles.ElementAt(a).Value[i][3]);
							}
							else
							{
								if (!epitopes_bindingcoreDic[alleles.ElementAt(a).Value[i][2]].Contains(alleles.ElementAt(a).Value[i][3]))
									epitopes_bindingcoreDic[alleles.ElementAt(a).Value[i][2]].Add(alleles.ElementAt(a).Value[i][3]);
							}
						}
						else if (version == "netMHCIIpan-3.0")
						{
							if (!epitopes_bindingcoreDic.ContainsKey(alleles.ElementAt(a).Value[i][2]))
							{
								epitopes_bindingcoreDic.Add(alleles.ElementAt(a).Value[i][2], new List<string>());
								epitopes_bindingcoreDic[alleles.ElementAt(a).Value[i][2]].Add(alleles.ElementAt(a).Value[i][5]);
							}
							else
							{
								if (!epitopes_bindingcoreDic[alleles.ElementAt(a).Value[i][2]].Contains(alleles.ElementAt(a).Value[i][5]))
									epitopes_bindingcoreDic[alleles.ElementAt(a).Value[i][2]].Add(alleles.ElementAt(a).Value[i][5]);
							}
						}

						else if (version == "iedb")
						{
							if (!epitopes_bindingcoreDic.ContainsKey(alleles.ElementAt(a).Value[i][5]))
							{
								epitopes_bindingcoreDic.Add(alleles.ElementAt(a).Value[i][5], new List<string>());
								epitopes_bindingcoreDic[alleles.ElementAt(a).Value[i][5]].Add(alleles.ElementAt(a).Value[i][4]);
							}
							else
							{
								if (!epitopes_bindingcoreDic[alleles.ElementAt(a).Value[i][5]].Contains(alleles.ElementAt(a).Value[i][4]))
									epitopes_bindingcoreDic[alleles.ElementAt(a).Value[i][5]].Add(alleles.ElementAt(a).Value[i][4]);
							}
						}
						else
							throw new Exception("0j03 version - " + version.ToString());
					}
				}
				foreach (KeyValuePair<string, List<string>> k in epitopes_bindingcoreDic)
				{
					tmp += k.Key + "(";
					foreach (string s in k.Value)
						tmp += s + ",";
					tmp = tmp.Substring(0, tmp.Length - 1);
					tmp += ")\t";
				}
				output.Add(tmp);
				//make the average row.
				List<double>[] averageList = new List<double>[1];
				double affinity = 0;

				foreach (KeyValuePair<string, List<string[]>> a in alleles)
				{
					averageList = new List<double>[a.Value.Count];

					for (int i = 0; i < averageList.Count(); i++)
						averageList[i] = new List<double>();
				}



				foreach (KeyValuePair<string, List<string[]>> a in alleles)
				{
					tmp = a.Key + "\t"; //display allele

					for (int i = 0; i < a.Value.Count; i++) //all epitope for 1 allele
					{

						//across
						if (version == "netMHCII-2.2")
						{
							tmp += a.Value[i][5] + "\t";
							if (!double.TryParse(a.Value[i][5], out affinity))
								throw new Exception("0j04 - " + tmp.ToString());
							else
								averageList[i].Add(affinity);
						}
						else if (version == "netMHCIIpan-3.0")
						{
							tmp += a.Value[i][7] + "\t";
							if (!double.TryParse(a.Value[i][7], out affinity))
								throw new Exception("0j04 - " + tmp.ToString());
							else
								averageList[i].Add(affinity);
						}
						else if (version == "iedb")
						{
							tmp += a.Value[i][6] + "\t";
							if (!double.TryParse(a.Value[i][6], out affinity))
								throw new Exception("0j04 - " + tmp.ToString());
							else
								averageList[i].Add(affinity);
						}
						else
							throw new Exception("0j03 version - " + version.ToString());
					}
					output.Add(tmp);

				}

				tmp = "average" + "\t";
				for (int i = 0; i < averageList.Count(); i++)
				{
					tmp += trimmedMean(averageList[i], 0.2).ToString() + "\t";
					//average[i] = int.Parse(l.Average ());
				}
				output.Add(tmp);
				output.Add("");
				l.write(output, path: (resultDir + "/" + version + "_results.tab"), delete: false);
				output.Clear();
				alleles.Clear();
				GC.Collect();
			}
		}

		public class parsemodified
		{
			public string seq { get; set;}
			public Dictionary<string, double> alleleAffinity { get; set;}
			public Dictionary<string,string>alleleCore { get; set; }
		}

		public List<parsemodified> parseModifiedResults()
		{
			Console.WriteLine("parseModifiedResults...");
			List<parsemodified> retArray = new List<parsemodified> ();
			string workingDir = root + "/netMHCII-2.2/modified_results";

			foreach (string n in Directory.GetDirectories(workingDir)) //this step loads everything into memory
			{					
				parsemodified p = new parsemodified ();
				Dictionary<string,double> affinity = new Dictionary<string,double> ();
				Dictionary<string,string> core = new Dictionary<string,string> ();
				string seq = "";

				foreach (string f in Directory.GetFiles(n))
				{
					List<string> input = l.read(path: f);
					string allele = f.Substring (f.LastIndexOf ("/") + 1).Replace (".result", ""); //get allele name
					//Console.WriteLine("Processing: " + protein + "/" + allele);

					for (int i = 1; i < input.Count; i++)
					{
						//filter.Add(input[i].Replace("\t", ","));
						if (input [i].Contains ("#")) {//ignore
						} else if (input [i].Contains ("---------")) {
						} else if (input [i].Contains (allele)) {
							//string temp = input [i].Replace ("\t", ",");
							string[] temp = input [i].Replace (" ", ",").Split (char.Parse (","));
							temp = temp.Where(w => w != "").ToArray();
							core.Add (allele, temp[3]);
							affinity.Add (allele, double.Parse(temp[5]));
							seq = (temp[2]);
							break;
						}
					}
				}

				p.alleleAffinity = affinity;
				p.seq = seq;
				p.alleleCore = core;
				retArray.Add (p);
				p = new parsemodified ();
				affinity = new Dictionary<string, double> ();
				core = new Dictionary<string, string> ();
				seq = "";
			}
			return retArray;
		}

		public double trimmedMean(List<double> numbers, double percent) //calcuatle trimmed mean of numbers, gotten from stack overflow
		{
			List<double> sorted = numbers;
			sorted.Sort();
			double sum = 0.0;
			double count = 0.0;
			int outlier = (int)Math.Floor(percent * numbers.Count() / 2);
			for (int i = outlier; i < sorted.Count() - outlier; i++)
			{
				sum += sorted[i];
				count++;
			}
			return (sum / count);
		}

		#region convert from tab to json
		/// <summary>
		/// takes the tab file derived from result files and summarize them into json files.
		/// </summary>
		public void parseTabToJson()
		{
			rawData rd = new rawData();
			List<Protein> proteinList = new List<Protein>();
			Protein proteins = new Protein();
			Epitope epitope = new Epitope();
			List<Epitope> epitopes = new List<Epitope>();
			List<Allele> allelesJson = new List<Allele>();
			Allele allele = new Allele();
			string proteinName, epitopeName, alleleName;
			double average;
			List<string> epitopeSeq = new List<string>();
			List<double> alleleAverage = new List<double>();
			Dictionary<string, List<double>> alleleIC50 = new Dictionary<string, List<double>>();
			List<string> temp = new List<string>();

			Dictionary<string, string> accessionToAllergen = new Dictionary<string, string>();
			foreach (string s in l.read(path: root + "/Protein_Accession_Dictionary.txt"))
			{
				if (!s.Contains(">") && !s.Contains("#") && s.Length < 2)
				{
					l.writeLog("[ERR]Improperly formated Protein_Accession_Dictionary.txt");
					throw new Exception("Improperly formated Protein_Accession_Dictionary.txt");
				}
				else if (s.Contains("#") || s.Length < 2)
				{
					//ignore comment and empty lines
				}
				else
					accessionToAllergen.Add(s.Substring(s.IndexOf(">") + 1), s.Substring(0, s.IndexOf(">")));
					//accessionToAllergen.Add(s.Substring(0, s.IndexOf(">")), s.Substring(s.IndexOf(">") + 1));
			}

			foreach (string file in Directory.GetFiles(resultDir))
			{
				if (file.Contains(".tab"))//result file
				{
					foreach (string input in l.read(path: file))
					{
						if (input.Contains("Antigen:") || input.Contains(".fasta"))//get protein name
						{
							if (!input.Contains(".fasta"))
								proteins.proteinName = input.Substring(input.IndexOf(": ") + 2);
							else
								proteins.proteinName = input.Substring(0, input.IndexOf("."));
							proteins.proteinAccession = accessionToAllergen[proteins.proteinName];//input.Substring(0, input.IndexOf("."));
							//remeber to change this for later.
							// string s = accessionToAllergen[proteins.proteinAccession];
							temp = l.read(path: root + "/Protein_Sequence/" + proteins.proteinAccession + ".fasta"); //!!!!!!!!!!!!!!!!!!!! change to .proteinaccession
							// proteins.proteinAccession = temp[0].Substring(l.NthOccurence(temp[0], Convert.ToChar("|"), 1) + 1, l.NthOccurence(temp[0], Convert.ToChar("|"), 2) - l.NthOccurence(temp[0], Convert.ToChar("|"), 1) - 1) ;
							temp.RemoveAt(0);
							proteins.proteinSeq = String.Join(String.Empty, temp.ToArray()).Replace(" ", "");
						}
						else if (input.Contains("#"))
						{
							//ignore comments and blank rows
						}
						else if (input == "" || input.Length <= 1)
						{
							//empty line to save everything
							//populate allelejson
							for (int i = 0; i < epitopeSeq.Count; i++)
							{
								epitope.epitopeSeq = epitopeSeq[i].Substring(0, epitopeSeq[i].IndexOf("("));
								epitope.average = alleleAverage[i];
								string epitopeList = epitopeSeq[i].Substring(epitopeSeq[i].IndexOf("(") + 1, epitopeSeq[i].IndexOf(")") - epitopeSeq[i].IndexOf("(") - 1);
								epitope.Binding_Core = epitopeList.Split(Convert.ToChar(",")).ToList();
								//now i need a List<allele>
								//allelesJson.Add()
								foreach (KeyValuePair<string, List<double>> k in alleleIC50)
								{
									allele = new Allele();
									allele.alleleName = k.Key;
									allele.IC50 = k.Value[i];
									allelesJson.Add(allele);
								}
								epitope.alleles = allelesJson;
								allelesJson = new List<Allele>();
								epitopes.Add(epitope);
								epitope = new Epitope();
							}
							proteins.epitopes = epitopes;
							epitopes = new List<Epitope>();
							proteinList.Add(proteins);
							proteins = new Protein();
							alleleIC50 = new Dictionary<string, List<double>>();
							epitopeSeq = new List<string>();
							alleleAverage = new List<double>();
						}
						else if (input.Substring(0, 1) == "\t")
						{
							//epitopes
							string tmp = input.Substring(1);
							do
							{
								if (tmp.Contains("\t"))
								{
									epitopeSeq.Add(tmp.Substring(0, tmp.IndexOf("\t")));
									tmp = tmp.Substring(tmp.IndexOf("\t") + 1);
								}
								else
								{
									break;
								}
							} while (true);
						}
						else if (input.Contains("average"))
						{
							string tmp = input.Substring(input.IndexOf("\t") + 1);
							do
							{
								if (tmp.Contains("\t"))
								{
									alleleAverage.Add(double.Parse(tmp.Substring(0, tmp.IndexOf("\t"))));
									tmp = tmp.Substring(tmp.IndexOf("\t") + 1);
								}
								else
								{
									break;
								}
							} while (true);
						}
						else if (input.Contains("HLA") || input.Contains("DRB") || input.Contains("result"))
						{
							//alleles
							//set allele name, ic50. 
							List<double> ic50 = new List<double>();
							string tmp = input;
							string n = input.Substring(0, input.IndexOf("."));
							tmp = input.Substring(input.IndexOf("\t") + 1);
							do
							{
								if (tmp.Contains("\t"))
								{
									ic50.Add(double.Parse(tmp.Substring(0, tmp.IndexOf("\t"))));
									tmp = tmp.Substring(tmp.IndexOf("\t") + 1);
								}
								else
									break;
							} while (true);
							alleleIC50.Add(n, ic50);
						}
						else
						{
							//???
						}
					}
					rd.all = proteinList;

					if (File.Exists(resultDir + "/" + file.Substring(file.LastIndexOf("/") + 1, file.LastIndexOf(".") - file.LastIndexOf("/") - 1) + ".json"))
						File.Move(resultDir + "/" + file.Substring(file.LastIndexOf("/") + 1, file.LastIndexOf(".") - file.LastIndexOf("/") - 1) + ".json", oldDirectory);
					using (StreamWriter f = File.AppendText(resultDir + "/" + file.Substring(file.LastIndexOf("/") + 1, file.LastIndexOf(".") - file.LastIndexOf("/") - 1) + ".json"))
					{
						JsonSerializer serializer = new JsonSerializer();
						serializer.Serialize(f, rd);
					}
					rd.all.Clear();
				}
			}
		}
		public class rawData
		{
			public List<Protein> all { get; set; }
		}
		public class Protein
		{
			public string proteinName { get; set; }
			public string proteinAccession { get; set; }
			public string proteinSeq { get; set; }
			public List<Epitope> epitopes { get; set; }
		}
		public class Epitope
		{
			public string epitopeSeq { get; set; }
			public List<string> Binding_Core { get; set; }
			public double average{ get; set; }
			public List<Population> populationCoverage { get; set; }
			public List<Allele> alleles { get; set; }
			public EpitopeParameters parameters { get; set; }
		}
		public class Allele
		{
			public string alleleName { get; set; }
			public double IC50 { get; set; }
		}
		public class Population
		{
			public string population { get; set; }
			public double coverage { get; set; }
			public double maxCoverage { get; set; }
			public bool madeCutoff { get; set; }
			public double Ranking_Score { get; set; }
			public bool Present_In_NetMHC2 { get; set; }
			public bool Present_In_NetMHC3 { get; set; }
			public bool Present_In_Tepitope { get; set; }
			public List<Allele> alleles { get; set; }
		}
		public class EpitopeParameters
		{
			public int Epitope_Start { get; set; }
			public int Epitope_Stop { get; set; }
			public double Molecular_Weight { get; set; }
			public int Extinction_Coefficient { get; set; }
			public double Isoelectric_point_pH { get; set; }
			public double Net_Charge_At_pH_7 { get; set; }
			public bool Innovagen_Solubility { get; set; }
			public double GRAVY_Score { get; set; }
			public ModifiedParameters Modified_Sequence_Parameters { get; set; }
		}
		public class ModifiedParameters
		{
			public string Modified_Sequence { get; set; }
			public bool Innovagen_Solubility { get; set; }
			public double GRAVY_Score { get; set; }
			public double Isoelectric_point_pH { get; set; }
			public double Net_Charge_At_pH_7 { get; set; }
			public List<Allele> Alleles{ get; set; }
			public List<string> Modified_Cores { get; set;}
			public double Modified_Average { get; set; }
		}

		#endregion
		#endregion

		/// <summary>
		/// function used to assign json to global scope.
		/// </summary>
		#region read json

		public void readJsonIntoMemory()
		{
			foreach (string s in Directory.GetFiles(resultDir)) //read in the 3 json files.
			{
				if (s.Contains(".json"))
				{
					if (s.Contains("netMHCII-2.2_results.json"))
					{
						using (StreamReader f = File.OpenText(s))
						{
							JsonSerializer serializer = new JsonSerializer();
							deserializedDatas[0] = (rawData)serializer.Deserialize(f, typeof(rawData));
						}
					}
					else if (s.Contains("netMHCIIpan-3.0_results.json"))
					{
						using (StreamReader f = File.OpenText(s))
						{
							JsonSerializer serializer = new JsonSerializer();
							deserializedDatas[1] = (rawData)serializer.Deserialize(f, typeof(rawData));
						}
					}
					else if (s.Contains("iedb_results.json"))
					{
						using (StreamReader f = File.OpenText(s))
						{
							JsonSerializer serializer = new JsonSerializer();
							deserializedDatas[2] = (rawData)serializer.Deserialize(f, typeof(rawData));
						}
					}
				}
			}

		}

		public void readAllDataIntoMemory()
		{
			Console.WriteLine ("Loading data into memory");
			foreach (string s in Directory.GetFiles(resultDir)) //read in the 3 json files.
			{
				if (s.Contains(".json"))
				{
					if (s.Contains("netMHCII-2.2_results_param.json"))
					{
						using (StreamReader f = File.OpenText(s))
						{
							JsonSerializer serializer = new JsonSerializer();
							deserializedDatas[0] = (rawData)serializer.Deserialize(f, typeof(rawData));
						}
					}
					else if (s.Contains("netMHCIIpan-3.0_results_param.json"))
					{
						using (StreamReader f = File.OpenText(s))
						{
							JsonSerializer serializer = new JsonSerializer();
							deserializedDatas[1] = (rawData)serializer.Deserialize(f, typeof(rawData));
						}
					}
					else if (s.Contains("iedb_results_param.json"))
					{
						using (StreamReader f = File.OpenText(s))
						{
							JsonSerializer serializer = new JsonSerializer();
							deserializedDatas[2] = (rawData)serializer.Deserialize(f, typeof(rawData));
						}
					}
				}
			}

		}
		#endregion

		/// <summary>
		/// function used to generate the heatmap and graphs for the prediction.
		/// </summary>
		#region generation Heatmap and Graph
		public void makeFigure()
		{
			#region generate heatmap
			Directory.CreateDirectory(resultDir + "/heatmaps");
			for (int i = 0; i < 2; i++)
			{
				Directory.CreateDirectory(resultDir + "/heatmaps/netMHCII_Ver." + (i+2).ToString());
				foreach (Protein protein in deserializedDatas[i].all)
				{
					using (Bitmap f = new Bitmap(protein.epitopes.Count, protein.epitopes[i].alleles.Count))
					{
						for (int column = 0; column < protein.epitopes.Count; column++)
							for (int row = 0; row < protein.epitopes[column].alleles.Count; row++)
								f.SetPixel(column, row, heatMapColorInc(protein.epitopes[column].alleles[row].IC50, 50, 500, 1000));

						f.Save(resultDir + "/heatmaps/netMHCII_Ver." + (i+2).ToString() + "/" + protein.proteinName + ".png", System.Drawing.Imaging.ImageFormat.Png);
					}
				}
			}

			Directory.CreateDirectory(resultDir + "/heatmaps/TEPITOPE");
			foreach (Protein protein in deserializedDatas[2].all)
			{
				using (Bitmap f = new Bitmap(protein.epitopes.Count, protein.epitopes[3].alleles.Count))
				{
					for (int column = 0; column < protein.epitopes.Count; column++)
						for (int row = 0; row < protein.epitopes[column].alleles.Count; row++)
							f.SetPixel(column, row, heatMapColorDsc(protein.epitopes[column].alleles[row].IC50, 1, 1.5, 2));

					f.Save(resultDir + "/heatmaps/TEPITOPE/" + protein.proteinName + ".png", System.Drawing.Imaging.ImageFormat.Png);
				}
			}
			#endregion

			#region generate graphs
			Directory.CreateDirectory(resultDir + "/charts");
			//foreach (Protein protein in deserializedDatas[0].all)
			for (int i = 0; i < deserializedDatas[0].all.Count; i++)
			{
				Dictionary<double, double> d = new Dictionary<double, double>();
				Dictionary<double, double> d2 = new Dictionary<double, double>();
				Dictionary<double, double> d3 = new Dictionary<double, double>();
				// Dictionary<double, double> d3 = new Dictionary<double, double>();

				for (int epitope = 0; epitope < deserializedDatas[0].all[i].epitopes.Count; epitope++)
				{ //pair average to key
					d.Add(double.Parse(epitope.ToString()), deserializedDatas[0].all[i].epitopes[epitope].average);
				}
				for (int epitope = 0; epitope < deserializedDatas[1].all[i].epitopes.Count; epitope++)
				{ //pair average to key
					d2.Add(double.Parse(epitope.ToString()), deserializedDatas[1].all[i].epitopes[epitope].average);
				}
				for (int epitope = 0; epitope < deserializedDatas[2].all[i].epitopes.Count; epitope++)
				{ //pair average to key
					d3.Add(double.Parse(epitope.ToString()), deserializedDatas[2].all[i].epitopes[epitope].average);//5000 * Math.Pow(2.718, (2.3 * deserializedDatas[2].all[i].epitopes[epitope].average)));
				}

				ChartUsingOxyPlot(d, d2,d3, deserializedDatas[0].all[i].proteinName);
				d.Clear();
				d2.Clear();
				d3.Clear();
				GC.Collect();
			}
			#endregion
		}
		/// <summary>
		/// return a color based on a number, min, middle and max
		/// </summary>
		/// <param name="value">a number</param>
		/// <param name="min">Minimum of number set</param>
		/// <param name="mid">Middle of number set</param>
		/// <param name="max">Maximum of number set</param>
		/// <returns></returns>
		public Color heatMapColorInc(double value, double min, double mid, double max)
		{
			//this is adapter from code obtained from stackoverflow available from http://stackoverflow.com/questions/5350235/how-to-generate-heat-maps-given-the-points
			//init
			int rOffset, gOffset, bOffset, deltaR, deltaG, deltaB, r, g, b;
			double val;
			//set the three color color spectrum starting colors
			Color firstColour = Color.Green;
			Color middleColour = Color.Yellow;
			Color secondColour = Color.Red;

			//for values less than the middle, color will be gradient from first color to middle color.
			if (value <= min)
				return firstColour;
			else if (value == mid) //middle number = middlecolor
				return middleColour;
			else if (value >= max)//largest number = second color
				return secondColour;
			else if (value < mid) //if value is smaller than mid = gradient between first and middle color
			{
				rOffset = Math.Max(firstColour.R, middleColour.R);
				gOffset = Math.Max(firstColour.G, middleColour.G);
				bOffset = Math.Max(firstColour.B, middleColour.B);

				deltaR = Math.Abs(middleColour.R - firstColour.R);
				deltaG = Math.Abs(middleColour.G - firstColour.G);
				deltaB = Math.Abs(middleColour.B - firstColour.B);

				val = (value - min) / (mid - min);
				r = rOffset - Convert.ToByte(deltaR * (1 - val));
				g = gOffset - Convert.ToByte(deltaG * (1 - val));
				b = bOffset - Convert.ToByte(deltaB * (1 - val));

				return Color.FromArgb(255, r, g, b);
			}
			else if (value > mid)//if value if greater than mid = gradient between middle and 2nd color
			{
				rOffset = Math.Max(middleColour.R, secondColour.R);
				gOffset = Math.Max(middleColour.G, secondColour.G);
				bOffset = Math.Max(middleColour.B, secondColour.B);

				deltaR = Math.Abs(middleColour.R - secondColour.R);
				deltaG = Math.Abs(middleColour.G - secondColour.G);
				deltaB = Math.Abs(middleColour.B - secondColour.B);

				val = (value - mid) / (max - mid);
				r = rOffset - Convert.ToByte(deltaR * (val));
				g = gOffset - Convert.ToByte(deltaG * (val));
				b = bOffset - Convert.ToByte(deltaB * (val));

				return Color.FromArgb(255, r, g, b);
			}
			else //shouldnt happen, error = white
				return Color.White;

		}
		/// <summary>
		/// return a color based on a number, min, middle and max
		/// </summary>
		/// <param name="value">a number</param>
		/// <param name="min">Minimum of number set</param>
		/// <param name="mid">Middle of number set</param>
		/// <param name="max">Maximum of number set</param>
		/// <returns></returns>
		public Color heatMapColorDsc(double value, double min, double mid, double max)
		{
			//this is adapter from code obtained from stackoverflow available from http://stackoverflow.com/questions/5350235/how-to-generate-heat-maps-given-the-points
			//init
			int rOffset, gOffset, bOffset, deltaR, deltaG, deltaB, r, g, b;
			double val;
			//set the three color color spectrum starting colors
			Color firstColour = Color.Red;
			Color middleColour = Color.Yellow;
			Color secondColour = Color.Green;

			//for values less than the middle, color will be gradient from first color to middle color.
			if (value <= min)
				return firstColour;
			else if (value == mid) //middle number = middlecolor
				return middleColour;
			else if (value >= max)//largest number = second color
				return secondColour;
			else if (value < mid) //if value is smaller than mid = gradient between first and middle color
			{
				rOffset = Math.Max(firstColour.R, middleColour.R);
				gOffset = Math.Max(firstColour.G, middleColour.G);
				bOffset = Math.Max(firstColour.B, middleColour.B);

				deltaR = Math.Abs(middleColour.R - firstColour.R);
				deltaG = Math.Abs(middleColour.G - firstColour.G);
				deltaB = Math.Abs(middleColour.B - firstColour.B);

				val = (value - min) / (mid - min);
				r = rOffset - Convert.ToByte(deltaR * (1 - val));
				g = gOffset - Convert.ToByte(deltaG * (1 - val));
				b = bOffset - Convert.ToByte(deltaB * (1 - val));

				return Color.FromArgb(255, r, g, b);
			}
			else if (value > mid)//if value if greater than mid = gradient between middle and 2nd color
			{
				rOffset = Math.Max(middleColour.R, secondColour.R);
				gOffset = Math.Max(middleColour.G, secondColour.G);
				bOffset = Math.Max(middleColour.B, secondColour.B);

				deltaR = Math.Abs(middleColour.R - secondColour.R);
				deltaG = Math.Abs(middleColour.G - secondColour.G);
				deltaB = Math.Abs(middleColour.B - secondColour.B);

				val = (value - mid) / (max - mid);
				r = rOffset - Convert.ToByte(deltaR * (val));
				g = gOffset - Convert.ToByte(deltaG * (val));
				b = bOffset - Convert.ToByte(deltaB * (val));

				return Color.FromArgb(255, r, g, b);
			}
			else //shouldnt happen, error = white
				return Color.White;

		}

		/// <summary>
		/// function used to call the charting form to draw and save the graph.
		/// </summary>
		/// <param name="value1">series 1 (2.2)</param>
		/// <param name="value2">series 2 (3.0)</param>
		/// <param name="value3">series 3 (tepitope)</param>
		/// <param name="seriesName"></param>
		public void charting(Dictionary<double, double> value1, Dictionary<double, double> value2, Dictionary<double,double> value3, string seriesName)
		{ 
			/*
			Chart ct = new Chart();

			System.Windows.Forms.DataVisualization.Charting.Series series1 = new System.Windows.Forms.DataVisualization.Charting.Series()
			{
				Name = "Average IC50, netmhc2.2",
				Color = Color.CornflowerBlue,
				BorderWidth = 3,
				IsVisibleInLegend = true,
				IsXValueIndexed = true,
				Enabled = true,
				ChartType = SeriesChartType.Line
			};
			System.Windows.Forms.DataVisualization.Charting.Series series2 = new System.Windows.Forms.DataVisualization.Charting.Series()
			{
				Name = "Average IC50, netmhc3.0",
				Color = Color.Green,
				BorderWidth = 3,
				IsVisibleInLegend = true,
				IsXValueIndexed = true,
				Enabled = true,
				ChartType = SeriesChartType.Line
			};
			System.Windows.Forms.DataVisualization.Charting.Series series3 = new System.Windows.Forms.DataVisualization.Charting.Series()
			{
				Name = "Average score, tepitope",
				Color = Color.Red,
				BorderWidth = 3,
				IsVisibleInLegend = true,
				IsXValueIndexed = true,
				Enabled = true,
				ChartType = SeriesChartType.Line
			};

			for (int i = 0; i < Math.Min(value1.Count, Math.Min(value2.Count, value3.Count)); i++)
			{
				series1.Points.AddXY(i + 1, value1[i]);
				series2.Points.AddXY(i + 1, value2[i]);
				series3.Points.AddXY(i + 1, value3[i]);
			}
			series1.Name = "netMHCII2.2";
			series2.Name = "netMHCII3.0";
			series3.Name = "Tepitope";

			Directory.CreateDirectory(resultDir + "/graphs");
			//c.charting(new Series[] { series1, series2, series3}, resultDir + "/graphs/" + seriesName + ".png");
			ct.Height = 365;
			ct.Width = 1000; 
			ct.Series.Clear();

			ct.ChartAreas.Add("Prediction Score");

			foreach (System.Windows.Forms.DataVisualization.Charting.Series se in new System.Windows.Forms.DataVisualization.Charting.Series[] { series1, series2, series3 })
			{
				ct.Series.Add(se);
				se.IsVisibleInLegend = true;
				ct.Legends.Add(se.Name);
			}

			foreach (Legend le in ct.Legends)
				le.Enabled = true;

			Axis xaxis = ct.ChartAreas[0].AxisX;
			xaxis.Interval = 30;
			xaxis.Title = "Amino Acid Number";
			xaxis.LineColor = Color.Gray;
			Axis yaxis = ct.ChartAreas[0].AxisY;
			yaxis.Interval = 1000;
			//yaxis.Maximum = 4000;
			yaxis.Title = "IC50";
			// yaxis.IsLogarithmic = true;
			yaxis.Maximum = 8000;
			// yaxis.LogarithmBase = 5;
			yaxis.LineColor = Color.Gray;

			ct.Invalidate(); //force draw before saving
			ct.SaveImage(resultDir + "/graphs/" + seriesName + ".png", ChartImageFormat.Png);
			/*using (Bitmap bmp = new Bitmap(ct.ClientRectangle.Width, ct.ClientRectangle.Height))
            {
                ct.DrawToBitmap(bmp, ct.ClientRectangle);
                bmp.Save(resultDir + "/graphs/" + seriesName + ".png", ImageFormat.Png);
            }*/


			//series1.Dispose();
			//series2.Dispose();
			//series3.Dispose();
			//ct.Dispose();
		}
		public void ChartUsingOxyPlot(Dictionary<double, double> value1, Dictionary<double, double> value2, Dictionary<double, double> value3, string seriesName)
		{
			var myModel = new PlotModel { Title = "Epitope Prediction Score" };

			LineSeries ls1 = new LineSeries();
			LineSeries ls2 = new LineSeries();
			LineSeries ls3 = new LineSeries();
			LineSeries indicator = new LineSeries();

			for (int i = 0; i < Math.Min(value1.Count, Math.Min(value2.Count, value3.Count)); i++)
			{
				ls1.Points.Add(new OxyPlot.DataPoint(i + 1, value1[i]));
				ls2.Points.Add(new OxyPlot.DataPoint(i + 1, value2[i]));
				ls3.Points.Add(new OxyPlot.DataPoint(i + 1, value3[i]));
				indicator.Points.Add(new OxyPlot.DataPoint(i + 1, 1000));
			}
			//legends
			ls1.Title = "netMHCII2.2";
			ls2.Title = "netMHCII3.0";
			ls3.Title = "Tepitope";

			ls1.MarkerType = MarkerType.Circle;
			ls1.MarkerFill= OxyColor.FromRgb(34, 145, 6);
			ls2.MarkerType = MarkerType.Circle;
			ls2.MarkerFill = OxyColor.FromRgb(6, 20, 145);
			ls3.MarkerType = MarkerType.Circle;
			ls3.MarkerFill = OxyColor.FromRgb(255, 0, 0);

			ls1.Color = OxyColor.FromRgb(34, 145, 6);
			ls2.Color = OxyColor.FromRgb(6, 20, 145);
			ls3.Color = OxyColor.FromRgb(255, 0, 0);
			indicator.Color = OxyColor.FromRgb(0,0,0);

			myModel.IsLegendVisible = true;
			myModel.LegendFontSize = 16;
			myModel.LegendPlacement = LegendPlacement.Outside;
			myModel.LegendPosition = LegendPosition.TopRight;
			//myModel.Axes.Add(new oxy LinearAxis(OxyPlot.Axes.AxisPosition.Left, 0, 8000));
			myModel.Axes.Add(new OxyPlot.Axes.LogarithmicAxis { Title = "Average NetMHCII Binding Affinity", Position = OxyPlot.Axes.AxisPosition.Left, StartPosition = 0, EndPosition = 1, Maximum = 25000, Minimum = 10,  });
			myModel.Axes.Add(new OxyPlot.Axes.LinearAxis { Title = "Average Tepitope Score", Position = OxyPlot.Axes.AxisPosition.Right, StartPosition = 1, EndPosition = 0, Maximum = 6, Minimum = -4 });

			ls3.YAxisKey = "Secondary";
			myModel.Axes[1].Key = "Secondary";

			ls1.YAxisKey = "PRIMARY";
			ls2.YAxisKey = "PRIMARY";

			myModel.Axes[0].Key = "PRIMARY";



			myModel.Series.Add(ls1);
			myModel.Series.Add(ls2);
			myModel.Series.Add(ls3);
			myModel.Series.Add(indicator);


			//myModel.DefaultXAxis.Title = "Amino Acid Start Position";
			//myModel.DefaultYAxis.Title = "Binding Score";

			using (var stream = File.Create(resultDir + "/charts/" + seriesName + ".svg"))
			{
				var exporter = new SvgExporter() { Width = 1000, Height = 500 };
				exporter.Export(myModel, stream);
			}
		}
		#endregion

		/// <summary>
		/// Functions used to generate the final report
		/// </summary>
		#region Filter and output the final results
		//this function writes populationCoverage Data to the json.
		public void getPopulationData()
		{
			double cutoff;

			if (DS.cutoff > 0 && DS.cutoff < 1)
				cutoff = DS.cutoff;
			else
				cutoff = 0.10;
			for (int i = 0; i < 2; i++)
			{
				foreach (Protein protein in deserializedDatas[i].all)
				{
					double sum = 0;
					double maxCoverage = 0;
					//List<string[]> StrongFilterArray = new List<string[]>();
					//List<string[]> WeakFilterArray = new List<string[]>();
					foreach (Epitope epitope in protein.epitopes)
					{
						List<Population> populations = new List<Population>();
						List<Allele> alleles = new List<Allele>();
						foreach (string s in Directory.GetFiles(root + "/Allele_Frequency"))
						{
							Population population = new Population();
							Dictionary<string, double> popFrequency = new Dictionary<string, double>();
							//load the population frequency data
							foreach (string data in l.read(path: s))
							{
								try
								{
									if (!data.Contains("#") && data.Length > 1)
										popFrequency.Add(data.Substring(0, data.IndexOf("(")).Replace("HLA", "").Replace("-","").Replace("*", "").Replace(":", "").Replace(" ", ""), double.Parse(data.Substring(data.LastIndexOf("(") + 1, data.LastIndexOf(")") - data.LastIndexOf("(") - 1)));
								}
								catch (Exception e)
								{
									l.writeLog(string.Format("[ERR]Error Reading Allele Frequency Data: {0} Error: {1}", data, e.ToString()));
								}
							}

							foreach (Allele allele in epitope.alleles)
							{
								if (allele.IC50 <= 500)
								{
									if (popFrequency.ContainsKey(allele.alleleName.Replace("HLA", "").Replace("-", "").Replace("_","")))
									{
										if (popFrequency[allele.alleleName.Replace("HLA", "").Replace("-", "").Replace("_", "")] != 0)
											sum += popFrequency[allele.alleleName.Replace("HLA", "").Replace("-", "").Replace("_", "")];
										else
											sum += 0.00001;
									}
									else
										sum += 0.00001;
									alleles.Add(allele);
								}

								if (popFrequency.ContainsKey(allele.alleleName.Replace("HLA", "").Replace("-", "").Replace("_", "")))
								{
									if (popFrequency[allele.alleleName.Replace("HLA", "").Replace("-", "").Replace("_", "")] != 0)
										maxCoverage += popFrequency[allele.alleleName.Replace("HLA", "").Replace("-", "").Replace("_", "")];
									else
										maxCoverage += 0.00001;
								}
								else
									maxCoverage += 0.00001;
							}
							if (sum > maxCoverage * (1-cutoff))
							{
								population.madeCutoff = true;
							}
							else
								population.madeCutoff = false;
							population.population = s.Substring(s.LastIndexOf("/") + 1, s.LastIndexOf(".") - s.LastIndexOf("/") - 1);
							population.coverage = sum;
							population.alleles = alleles;
							alleles = new List<Allele>();
							population.maxCoverage = maxCoverage;
							populations.Add(population);
							sum = 0;
							maxCoverage = 0;
						}


						epitope.populationCoverage = populations;
						populations = new List<Population>();
					}
				}

				string tmpDir;
				if (i == 0)
					tmpDir = resultDir + "/netMHCII-2.2_results_param.json";
				else
					tmpDir = resultDir + "/netMHCIIpan-3.0_results_param.json";
				File.Delete(tmpDir);
				using (StreamWriter f = File.AppendText(tmpDir))
				{
					foreach (Protein p in deserializedDatas[i].all)
						foreach (Epitope e in p.epitopes)
							e.alleles = null;
					JsonSerializer serializer = new JsonSerializer();
					serializer.Serialize(f, deserializedDatas[i]);
				}
			}
			//tepitope
			foreach (Protein protein in deserializedDatas[2].all)
			{
				double sum = 0;
				double maxCoverage = 0;
				List<string[]> StrongFilterArray = new List<string[]>();
				List<string[]> WeakFilterArray = new List<string[]>();
				foreach (Epitope epitope in protein.epitopes)
				{
					List<Population> populations = new List<Population>();
					List<Allele> alleles = new List<Allele>();
					foreach (string s in Directory.GetFiles(root + "/Allele_Frequency"))
					{
						Population population = new Population();
						Dictionary<string, double> popFrequency = new Dictionary<string, double>();
						//load the population frequency data
						foreach (string data in l.read(path: s))
						{
							try
							{
								if (!data.Contains("#") && data.Length > 1)
									popFrequency.Add(data.Substring(0, data.IndexOf("(")).Replace("HLA", "").Replace("-", "").Replace("*", "").Replace(":", "").Replace(" ", ""), double.Parse(data.Substring(data.LastIndexOf("(") + 1, data.LastIndexOf(")") - data.LastIndexOf("(") - 1)));
							}
							catch (Exception e)
							{
								l.writeLog(string.Format("[ERR]Error Reading Allele Frequency Data: {0} Error: {1}", data, e.ToString()));
							}
						}

						foreach (Allele allele in epitope.alleles)
						{
							if (allele.IC50 >= 1)
							{
								if (popFrequency.ContainsKey(allele.alleleName.Replace("HLA", "").Replace("-", "").Replace("*", "").Replace(":", "")))
								{
									if (popFrequency[allele.alleleName.Replace("HLA", "").Replace("-", "").Replace("*", "").Replace(":", "")] != 0)
										sum += popFrequency[allele.alleleName.Replace("HLA", "").Replace("-", "").Replace("*", "").Replace(":", "")];
									else
										sum += 0.00001;
								}
								else
									sum += 0.00001;
								alleles.Add(allele);
							}

							if (popFrequency.ContainsKey(allele.alleleName.Replace("HLA", "").Replace("-", "").Replace("*", "").Replace(":", "")))
							{
								if (popFrequency[allele.alleleName.Replace("HLA", "").Replace("-", "").Replace("*", "").Replace(":", "")] != 0)
									maxCoverage += popFrequency[allele.alleleName.Replace("HLA", "").Replace("-", "").Replace("*", "").Replace(":", "")];
								else
									maxCoverage += 0.00001;
							}
							else
								maxCoverage += 0.00001;
						}
						if (sum > maxCoverage * (1 - cutoff))
						{
							population.madeCutoff = true;
						}
						else
							population.madeCutoff = false;
						population.population = s.Substring(s.LastIndexOf("/") + 1, s.LastIndexOf(".") - s.LastIndexOf("/") - 1);
						population.coverage = sum;
						population.alleles = alleles;
						alleles = new List<Allele>();
						population.maxCoverage = maxCoverage;
						//population.Final_Score = sum * epitope.average;
						populations.Add(population);
						sum = 0;
						maxCoverage = 0;
					}


					epitope.populationCoverage = populations;
					populations = new List<Population>();
				}
			}

			string td = resultDir + "/iedb_results_param.json";
			File.Delete(td);
			using (StreamWriter f = File.AppendText(td))
			{
				//foreach (Protein p in deserializedDatas[2].all)
				//    foreach (Epitope e in p.epitopes)
				//      e.alleles = null;
				JsonSerializer serializer = new JsonSerializer();
				serializer.Serialize(f, deserializedDatas[2]);
			}
		}
		//this function writes epitope parameters Data to the json.
		public void getEpitopeParameter()
		{
			//check if made cut off is true, then assign parameters
			for (int pr = 0; pr < deserializedDatas[0].all.Count; pr++) //foreach protein
			{
				for (int e = 0; e < deserializedDatas[0].all[pr].epitopes.Count; e++) //foreach epitope
				{
					if (deserializedDatas[0].all[pr].epitopes[e].epitopeSeq != deserializedDatas[1].all[pr].epitopes[e].epitopeSeq || deserializedDatas[0].all[pr].epitopes[e].epitopeSeq != deserializedDatas[2].all[pr].epitopes[e].epitopeSeq || deserializedDatas[1].all[pr].epitopes[e].epitopeSeq != deserializedDatas[2].all[pr].epitopes[e].epitopeSeq)
					{
						throw new Exception("0j919");
					}

					//make parameters

					for (int po = 0; po < deserializedDatas[0].all[pr].epitopes[e].populationCoverage.Count; po++) //foreach population
					{
						//generate a list for each population in json.
						if (deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].madeCutoff)
						{
							deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC2 = true;
							deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC2 = true;
							deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC2 = true;
						}
						else
						{
							deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC2 = false;
							deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC2 = false;
							deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC2 = false;
						}

						if (deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].madeCutoff)
						{
							deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC3 = true;
							deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC3 = true;
							deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC3 = true;
						}
						else
						{
							deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC3 = false;
							deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC3 = false;
							deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC3 = false;
						}

						if (deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].madeCutoff)
						{
							deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].Present_In_Tepitope = true;
							deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].Present_In_Tepitope = true;
							deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Present_In_Tepitope = true;
						}
						else
						{
							deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].Present_In_Tepitope = false;
							deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].Present_In_Tepitope = false;
							deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Present_In_Tepitope = false;
						}

						//calculate final score
						deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].Ranking_Score =
							((deserializedDatas[0].all[pr].epitopes[e].average * (1-(deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].coverage / deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].maxCoverage)) +
								deserializedDatas[1].all[pr].epitopes[e].average * (1-(deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].coverage / deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].maxCoverage))) / 2) /
								(Math.Abs(deserializedDatas[2].all[pr].epitopes[e].average + 4) * (deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].coverage / deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].maxCoverage));
						deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].Ranking_Score = deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].Ranking_Score;
						deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Ranking_Score = deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].Ranking_Score;
						//get parameters for those made cutoff. 

						//have to change this so that parameter isnt gotten for every single thing
						if (deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Present_In_Tepitope || deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC3 || deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC2)
						{
							if (deserializedDatas[0].all[pr].epitopes[e].parameters == null || deserializedDatas[1].all[pr].epitopes[e].parameters == null || deserializedDatas[2].all[pr].epitopes[e].parameters == null)
							{
								deserializedDatas[0].all[pr].epitopes[e].parameters = getInnovagenSolubility(deserializedDatas[0].all[pr].epitopes[e].epitopeSeq);
								deserializedDatas[0].all[pr].epitopes[e].parameters.Epitope_Start = e + 1;// deserializedDatas[0].all[pr].proteinSeq.IndexOf(deserializedDatas[0].all[pr].epitopes[e].epitopeSeq);
								deserializedDatas[0].all[pr].epitopes[e].parameters.Epitope_Stop = e + 15;// deserializedDatas[0].all[pr].epitopes[e].parameters.Epitope_Start + 15;
								if (!deserializedDatas[0].all[pr].epitopes[e].parameters.Innovagen_Solubility)
									deserializedDatas[0].all[pr].epitopes[e].parameters.Modified_Sequence_Parameters = optimize(deserializedDatas[0].all[pr].epitopes[e].epitopeSeq, deserializedDatas[0].all[pr].proteinSeq, deserializedDatas[0].all[pr].epitopes[e].parameters.Epitope_Start, deserializedDatas[0].all[pr].epitopes[e].parameters.Epitope_Stop);
								deserializedDatas[1].all[pr].epitopes[e].parameters = deserializedDatas[0].all[pr].epitopes[e].parameters;
								deserializedDatas[1].all[pr].epitopes[e].parameters.Epitope_Start = deserializedDatas[0].all[pr].epitopes[e].parameters.Epitope_Start;
								deserializedDatas[1].all[pr].epitopes[e].parameters.Epitope_Stop = deserializedDatas[0].all[pr].epitopes[e].parameters.Epitope_Stop;
								deserializedDatas[1].all[pr].epitopes[e].parameters.Modified_Sequence_Parameters = deserializedDatas[0].all[pr].epitopes[e].parameters.Modified_Sequence_Parameters;
								deserializedDatas[2].all[pr].epitopes[e].parameters = deserializedDatas[0].all[pr].epitopes[e].parameters;
								deserializedDatas[2].all[pr].epitopes[e].parameters.Epitope_Start = deserializedDatas[0].all[pr].epitopes[e].parameters.Epitope_Start;
								deserializedDatas[2].all[pr].epitopes[e].parameters.Epitope_Stop = deserializedDatas[0].all[pr].epitopes[e].parameters.Epitope_Stop;
								deserializedDatas[2].all[pr].epitopes[e].parameters.Modified_Sequence_Parameters = deserializedDatas[0].all[pr].epitopes[e].parameters.Modified_Sequence_Parameters;
							}
						}
					}

				}
			}

			File.Delete(resultDir + "/netMHCII-2.2_results_param.json");
			using (StreamWriter f = File.AppendText(resultDir + "/netMHCII-2.2_results_param.json"))
			{
				JsonSerializer serializer = new JsonSerializer();
				serializer.Serialize(f, deserializedDatas[0]);
			}
			Thread.Sleep(1000);
			File.Delete(resultDir + "/netMHCIIpan-3.0_results_param.json");
			using (StreamWriter f = File.AppendText(resultDir + "/netMHCIIpan-3.0_results_param.json"))
			{
				JsonSerializer serializer = new JsonSerializer();
				serializer.Serialize(f, deserializedDatas[1]);
			}
			Thread.Sleep(1000);
			File.Delete(resultDir + "/iedb_results_param.json");
			using (StreamWriter f = File.AppendText(resultDir + "/iedb_results_param.json"))
			{
				JsonSerializer serializer = new JsonSerializer();
				serializer.Serialize(f, deserializedDatas[2]);
			}

		}
		//output the results into tab
		public void outputResult()
		{
			Directory.CreateDirectory(resultDir + "/Result_Summary");
			List<string>[] outputs = new List<string>[deserializedDatas[0].all[0].epitopes[0].populationCoverage.Count];
			//List<string>[] condense = new List<string>[deserializedDatas[0].all[0].epitopes[0].populationCoverage.Count];
			//List<finalJson>[] result = new List<finalJson>[deserializedDatas[0].all[0].epitopes[0].populationCoverage.Count];
			List<Epitope> ep; List<Population> pop;
			int[] lastStart = new int[outputs.Length];
			int[] lastPosition = new int[outputs.Length];
			string[] lastcore0 = new string[outputs.Length];
			string[] lastcore1 = new string[outputs.Length];
			string[] lastcore2 = new string[outputs.Length];
			string[] combinedEpitope = new string[outputs.Length];
			bool[] lastcoverage0 = new bool[outputs.Length];
			bool[] lastcoverage1 = new bool[outputs.Length];
			bool[] lastcoverage2 = new bool[outputs.Length];
			// bool[] print = new bool[outputs.Length];
			Epitope[] epitope = new Epitope[outputs.Length];
			Epitope[,] lastEpitope = new Epitope[outputs.Length, 3];

			for (int i = 0; i < outputs.Length; i++)
			{
				lastStart[i] = -99;
				lastEpitope[i, 0] = new Epitope();
				lastEpitope[i, 0].parameters = new EpitopeParameters();
				lastEpitope[i, 1] = new Epitope();
				lastEpitope[i, 1].parameters = new EpitopeParameters();
				lastEpitope[i, 2] = new Epitope();
				lastEpitope[i, 2].parameters = new EpitopeParameters();
				epitope[i] = new Epitope();
				epitope[i].parameters = new EpitopeParameters();
				//  print[i] = false;
				outputs[i] = new List<string>();
				outputs[i].Add("#" + deserializedDatas[0].all[0].epitopes[0].populationCoverage[i].population);
				//condense[i] = new List<string>();
				//condense[i].Add("#" + deserializedDatas[0].all[0].epitopes[0].populationCoverage[i].population);
				//condense[i] = new List<string>();
			}

			for (int pr = 0; pr < deserializedDatas[0].all.Count; pr++) //foreach protein
			{
				for (int i = 0; i < outputs.Length; i++)
				{
					outputs[i].Add("Antigen: " + deserializedDatas[0].all[pr].proteinName);
					//bla g 1.1, 	coverage	binding core	Present in NetMHC2.2	Present in NetMHCpan3.0	Present in Tepitope	Start	Stop	Molecular Weight	Extinction Coefficient	Iso-electric point	Net charge at pH 7	Innovagen Solubility	GRAVY score	Soluble Modification
					outputs[i].Add("Epitope Sequence\tOverall Score\tBinding Core(NetMHC2, NetMHC3, Tepitope)\tPredicted Score(NetMHC2, NetMHC3, Tepitope)\tPopulation Coverage(NetMHC2/Max Coverage, NetMHC3/Max Coverage, Tepitope/Max Coverage)\tPredicted by(NetMHC2, NetMHC3, Tepitope)\tStart\tStop\tMolecular Weight\tExtinction Coefficient\tIso-electric point\tNet Charge at pH7\tSolubility\tGRAVY Score\tModified Sequence\tModified GRAVY Score\tModified Solubility\tModified NetMHC2.2 Score)");

					//condense[i].Add("Antigen: " + deserializedDatas[0].all[pr].proteinName);
					//bla g 1.1, 	coverage	binding core	Present in NetMHC2.2	Present in NetMHCpan3.0	Present in Tepitope	Start	Stop	Molecular Weight	Extinction Coefficient	Iso-electric point	Net charge at pH 7	Innovagen Solubility	GRAVY score	Soluble Modification
					//condense[i].Add("Epitope Sequence\tOverall Score\tBinding Core\tStart\tStop\tMolecular Weight\tExtinction Coefficient\tIso-electric point\tNet Charge at pH7\tSolubility\tGRAVY Score\tModified Sequence\tModified GRAVY Score\t Modified Solubility)");
					// print[i] = true;
				}
				ep = deserializedDatas[0].all[pr].epitopes;
				for (int e = 0; e < ep.Count; e++) //foreach epitope
				{
					pop = deserializedDatas[0].all[pr].epitopes[e].populationCoverage;
					for (int po = 0; po < pop.Count; po++) //foreach population
					{
						if (pop[po].Present_In_NetMHC2 || pop[po].Present_In_NetMHC3 || pop[po].Present_In_Tepitope)
						{
							//got everything
							if (ep[e].parameters.Modified_Sequence_Parameters != null)
							{
								//outputs[i].Add("\t\tStart\tStop\tMolecular Weight\tExtinction Coefficient\tIso-electric point\tNet Charge at pH7\tSolubility\tGRAVY Score\tModified Sequence(If Necessary)(GRAVY Score");
								outputs[po].Add(String.Format("{0}\t{1}\t({2}|{3}|{4})\t({5}|{6}|{7})\t({8}/{11}|{9}/{12}|{10}/{13})\t({14}|{15}|{16})\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}\t{25}\t{26}\t{27}\t{28}\t{29}",
								                              ep[e].epitopeSeq, String.Format("{0:0.00}", deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].Ranking_Score),
								                              string.Join(",", deserializedDatas[0].all[pr].epitopes[e].Binding_Core), string.Join(",", deserializedDatas[1].all[pr].epitopes[e].Binding_Core), string.Join(",", deserializedDatas[2].all[pr].epitopes[e].Binding_Core),
								                              string.Format("{0:0.00}", deserializedDatas[0].all[pr].epitopes[e].average), string.Format("{0:0.00}", deserializedDatas[1].all[pr].epitopes[e].average), string.Format("{0:0.00}", deserializedDatas[2].all[pr].epitopes[e].average),
								                              string.Format("{0:0.00000}", deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].coverage), string.Format("{0:0.00000}", deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].coverage), string.Format("{0:0.00000}", deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].coverage),
								                              string.Format("{0:0.00}", deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].maxCoverage), string.Format("{0:0.00}", deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].maxCoverage), string.Format("{0:0.00}", deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].maxCoverage),
								                              pop[po].Present_In_NetMHC2.ToString().Replace("True", "Yes").Replace("False", "No"), pop[po].Present_In_NetMHC3.ToString().Replace("True", "Yes").Replace("False", "No"), pop[po].Present_In_Tepitope.ToString().Replace("True", "Yes").Replace("False", "No"),
								                              ep[e].parameters.Epitope_Start, ep[e].parameters.Epitope_Stop, ep[e].parameters.Molecular_Weight, ep[e].parameters.Extinction_Coefficient, ep[e].parameters.Isoelectric_point_pH,
									ep[e].parameters.Net_Charge_At_pH_7, ep[e].parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), ep[e].parameters.GRAVY_Score, ep[e].parameters.Modified_Sequence_Parameters.Modified_Sequence, ep[e].parameters.Modified_Sequence_Parameters.GRAVY_Score, ep[e].parameters.Modified_Sequence_Parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), ep[e].parameters.Modified_Sequence_Parameters.Modified_Average, string.Join(",",ep[e].parameters.Modified_Sequence_Parameters.Modified_Cores)
								                              ));
							}
							else
							{
								outputs[po].Add(String.Format("{0}\t{1}\t({2}|{3}|{4})\t({5}|{6}|{7})\t({8}/{11}|{9}/{12}|{10}/{13})\t({14}|{15}|{16})\t{17}\t{18}\t{19}\t{20}\t{21}\t{22}\t{23}\t{24}",
								                              ep[e].epitopeSeq, String.Format("{0:0.00}", deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].Ranking_Score),
								                              string.Join(",", deserializedDatas[0].all[pr].epitopes[e].Binding_Core), string.Join(",", deserializedDatas[1].all[pr].epitopes[e].Binding_Core), string.Join(",", deserializedDatas[2].all[pr].epitopes[e].Binding_Core),
								                              string.Format("{0:0.00}", deserializedDatas[0].all[pr].epitopes[e].average), string.Format("{0:0.00}", deserializedDatas[1].all[pr].epitopes[e].average), string.Format("{0:0.00}", deserializedDatas[2].all[pr].epitopes[e].average),
								                              string.Format("{0:0.00000}", deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].coverage), string.Format("{0:0.00000}", deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].coverage), string.Format("{0:0.00000}", deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].coverage),
								                              string.Format("{0:0.00}", deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].maxCoverage), string.Format("{0:0.00}", deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].maxCoverage), string.Format("{0:0.00}", deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].maxCoverage),
								                              pop[po].Present_In_NetMHC2.ToString().Replace("True", "Yes").Replace("False", "No"), pop[po].Present_In_NetMHC3.ToString().Replace("True", "Yes").Replace("False", "No"), pop[po].Present_In_Tepitope.ToString().Replace("True", "Yes").Replace("False", "No"),
								                              ep[e].parameters.Epitope_Start, ep[e].parameters.Epitope_Stop, ep[e].parameters.Molecular_Weight, ep[e].parameters.Extinction_Coefficient, ep[e].parameters.Isoelectric_point_pH,
								                              ep[e].parameters.Net_Charge_At_pH_7, ep[e].parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), ep[e].parameters.GRAVY_Score
								                              ));
							}

							//no this doesnt work.
							//try to make a variable to stored the epitopes in a list. 
							int temp = e - 1;
							if (temp < 0) temp = 0;
							//filtered down
							/*if (
                                deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC2 == lastEpitope[po,0].populationCoverage[po].Present_In_NetMHC2 &&
                                deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC3 == lastEpitope[po, 1].populationCoverage[po].Present_In_NetMHC3 &&
                                deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Present_In_Tepitope == lastEpitope[po, 2].populationCoverage[po].Present_In_Tepitope &&)
                            {*/
							if (deserializedDatas[2].all[pr].epitopes[e].parameters.Epitope_Start < 1) //beginning
							{
								epitope[po] = new Epitope();
								epitope[po].epitopeSeq = deserializedDatas[2].all[pr].epitopes[e].epitopeSeq;
								//epitope[po].Binding_Core = deserializedDatas[1].all[pr].epitopes[e].Binding_Core;
								epitope[po].average = deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Ranking_Score;
								epitope[po].parameters = new EpitopeParameters();
								epitope[po].parameters.Epitope_Stop = deserializedDatas[2].all[pr].epitopes[e].parameters.Epitope_Stop;
								// print[po] = false;
								lastEpitope[po, 0] = deserializedDatas[0].all[pr].epitopes[e];
								lastEpitope[po, 1] = deserializedDatas[1].all[pr].epitopes[e];
								lastEpitope[po, 2] = deserializedDatas[2].all[pr].epitopes[e];
								lastStart[po] = deserializedDatas[2].all[pr].epitopes[e].parameters.Epitope_Start;


								if (ep[e].parameters.Modified_Sequence_Parameters != null)
								{
									//outputs[i].Add("\t\tStart\tStop\tMolecular Weight\tExtinction Coefficient\tIso-electric point\tNet Charge at pH7\tSolubility\tGRAVY Score\tModified Sequence(If Necessary)(GRAVY Score");
									/*condense[po].Add(String.Format("{0}\t{1}\t{2}{3}{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}",
									                               ep[e].epitopeSeq, String.Format("{0:0.00}", ep[e].populationCoverage[po].Ranking_Score),
									                               "", string.Join(",", deserializedDatas[1].all[pr].epitopes[e].Binding_Core), "",
									                               ep[e].parameters.Epitope_Start, ep[e].parameters.Epitope_Stop, ep[e].parameters.Molecular_Weight, ep[e].parameters.Extinction_Coefficient, ep[e].parameters.Isoelectric_point_pH,
									                               ep[e].parameters.Net_Charge_At_pH_7, ep[e].parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), ep[e].parameters.GRAVY_Score, ep[e].parameters.Modified_Sequence_Parameters.Modified_Sequence, ep[e].parameters.Modified_Sequence_Parameters.GRAVY_Score, ep[e].parameters.Modified_Sequence_Parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor")
									                               ));*/
								}
								else
								{
									/*condense[po].Add(String.Format("{0}\t{1}\t{2}{3}{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}",
									                               ep[e].epitopeSeq, String.Format("{0:0.00}", ep[e].populationCoverage[po].Ranking_Score),
									                               "", string.Join(",", deserializedDatas[1].all[pr].epitopes[e].Binding_Core), "",
									                               ep[e].parameters.Epitope_Start, ep[e].parameters.Epitope_Stop, ep[e].parameters.Molecular_Weight, ep[e].parameters.Extinction_Coefficient, ep[e].parameters.Isoelectric_point_pH,
									                               ep[e].parameters.Net_Charge_At_pH_7, ep[e].parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), ep[e].parameters.GRAVY_Score
									                               ));*/
								}
							}
							else if (deserializedDatas[2].all[pr].epitopes[e].parameters.Epitope_Start == lastEpitope[po, 2].parameters.Epitope_Start + 1)
							{
								//if (print[po] == false)//deserializedDatas[2].all[pr].epitopes[e].parameters.Epitope_Start == 1)
								// if (!print[po])
								//if (deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Ranking_Score < lastEpitope[po, 2].populationCoverage[po].Ranking_Score)
								// {
								//condense[po].RemoveAt(condense[po].Count - 1);//hax

								if (epitope[po].epitopeSeq == null)
									epitope[po].epitopeSeq = deserializedDatas[0].all[pr].epitopes[e].epitopeSeq;
								else
									epitope[po].epitopeSeq = epitope[po].epitopeSeq + deserializedDatas[0].all[pr].epitopes[e].epitopeSeq.Substring(14);
								/*
                                if (epitope[po].Binding_Core == null)// string.Format("({0}, {1}, {2})", deserializedDatas[0].all[pr].epitopes[e].Binding_Core, deserializedDatas[1].all[pr].epitopes[e].Binding_Core, deserializedDatas[2].all[pr].epitopes[e].Binding_Core)))
                                {
                                    epitope[po].Binding_Core += deserializedDatas[1].all[pr].epitopes[e].Binding_Core + ",";//string.Format("({0}, {1}, {2})", deserializedDatas[0].all[pr].epitopes[e].Binding_Core, deserializedDatas[1].all[pr].epitopes[e].Binding_Core, deserializedDatas[2].all[pr].epitopes[e].Binding_Core) + ",";
                                }
                                else if (!epitope[po].Binding_Core.Contains(deserializedDatas[1].all[pr].epitopes[e].Binding_Core))
                                {
                                    epitope[po].Binding_Core += deserializedDatas[1].all[pr].epitopes[e].Binding_Core + ",";
                                }*/
								epitope[po].Binding_Core = deserializedDatas[1].all[pr].epitopes[e].Binding_Core;
								epitope[po].average = (epitope[po].average + deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].Ranking_Score) / 2;
								// print[po] = true;
								lastEpitope[po, 0] = deserializedDatas[0].all[pr].epitopes[e];
								lastEpitope[po, 1] = deserializedDatas[1].all[pr].epitopes[e];
								lastEpitope[po, 2] = deserializedDatas[2].all[pr].epitopes[e];

								// epitope[po].parameters = getInnovagenSolubility(epitope[po].epitopeSeq);
								//epitope[po].parameters.Epitope_Start = lastStart[po];
								// epitope[po].parameters.Epitope_Stop = lastEpitope[po, 0].parameters.Epitope_Stop;
								//  if (!epitope[po].parameters.Innovagen_Solubility)
								//      epitope[po].parameters.Modified_Sequence_Parameters = optimize(epitope[po].epitopeSeq, "KKKK", 99, 94);
								//if (epitope[po].parameters.Modified_Sequence_Parameters != null) //take out 5-13.
								//{
								//outputs[i].Add("\t\tStart\tStop\tMolecular Weight\tExtinction Coefficient\tIso-electric point\tNet Charge at pH7\tSolubility\tGRAVY Score\tModified Sequence(If Necessary)(GRAVY Score");
								/*condense[po].Add(String.Format("{0}\t{1}\t{2}\t{3}\t{4}",
								                               epitope[po].epitopeSeq, String.Format("{0:0.00}", epitope[po].average),
								                               string.Join(",",epitope[po].Binding_Core), lastStart[po], lastStart[po] + epitope[po].epitopeSeq.Length - 1));//, epitope[po].parameters.Molecular_Weight, epitope[po].parameters.Extinction_Coefficient, epitope[po].parameters.Isoelectric_point_pH,
								// print[po] = true;                                                                                                  //epitope[po].parameters.Net_Charge_At_pH_7, epitope[po].parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), epitope[po].parameters.GRAVY_Score, epitope[po].parameters.Modified_Sequence_Parameters.Modified_Sequence, epitope[po].parameters.Modified_Sequence_Parameters.GRAVY_Score, epitope[po].parameters.Modified_Sequence_Parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor")
								// ));
								/*}
                                                                                                                                 else
                                                                                                                                 {
                                                                                                                                     condense[po].Add(String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}",
                                                                                                                                         epitope[po].epitopeSeq, String.Format("{0:0.00}", epitope[po].average),
                                                                                                                                         epitope[po].Binding_Core,
                                                                                                                                         epitope[po].parameters.Epitope_Start, epitope[po].parameters.Epitope_Stop, epitope[po].parameters.Molecular_Weight, epitope[po].parameters.Extinction_Coefficient, epitope[po].parameters.Isoelectric_point_pH,
                                                                                                                                         epitope[po].parameters.Net_Charge_At_pH_7, epitope[po].parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), epitope[po].parameters.GRAVY_Score
                                                                                                                                         ));
                                                                                                                                 }*/
								//}
							}
							else
							{
								if (ep[e].parameters.Modified_Sequence_Parameters != null)
								{
									//outputs[i].Add("\t\tStart\tStop\tMolecular Weight\tExtinction Coefficient\tIso-electric point\tNet Charge at pH7\tSolubility\tGRAVY Score\tModified Sequence(If Necessary)(GRAVY Score");
									/*condense[po].Add(String.Format("{0}\t{1}\t{2}{3}{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}",
									                               ep[e].epitopeSeq, String.Format("{0:0.00}", ep[e].populationCoverage[po].Ranking_Score),
									                               "", string.Join(",", deserializedDatas[1].all[pr].epitopes[e].Binding_Core), "",
									                               ep[e].parameters.Epitope_Start, ep[e].parameters.Epitope_Stop, ep[e].parameters.Molecular_Weight, ep[e].parameters.Extinction_Coefficient, ep[e].parameters.Isoelectric_point_pH,
									                               ep[e].parameters.Net_Charge_At_pH_7, ep[e].parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), ep[e].parameters.GRAVY_Score, ep[e].parameters.Modified_Sequence_Parameters.Modified_Sequence, ep[e].parameters.Modified_Sequence_Parameters.GRAVY_Score, ep[e].parameters.Modified_Sequence_Parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor")
									                               ));*/
								}
								else
								{
									/*condense[po].Add(String.Format("{0}\t{1}\t{2}{3}{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}",
									                               ep[e].epitopeSeq, String.Format("{0:0.00}", ep[e].populationCoverage[po].Ranking_Score),
									                               "", string.Join(",", deserializedDatas[1].all[pr].epitopes[e].Binding_Core), "",
									                               ep[e].parameters.Epitope_Start, ep[e].parameters.Epitope_Stop, ep[e].parameters.Molecular_Weight, ep[e].parameters.Extinction_Coefficient, ep[e].parameters.Isoelectric_point_pH,
									                               ep[e].parameters.Net_Charge_At_pH_7, ep[e].parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), ep[e].parameters.GRAVY_Score
									                               ));*/
								}
								// print[po] = true;
								epitope[po] = new Epitope();
								epitope[po].epitopeSeq = deserializedDatas[2].all[pr].epitopes[e].epitopeSeq;
								epitope[po].average = deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Ranking_Score;
								//epitope[po].Binding_Core = deserializedDatas[1].all[pr].epitopes[e].Binding_Core;
								epitope[po].parameters = new EpitopeParameters();
								epitope[po].parameters.Epitope_Stop = deserializedDatas[2].all[pr].epitopes[e].parameters.Epitope_Stop;
								//    print[po] = false;
								lastEpitope[po, 0] = deserializedDatas[0].all[pr].epitopes[e];
								lastEpitope[po, 1] = deserializedDatas[1].all[pr].epitopes[e];
								lastEpitope[po, 2] = deserializedDatas[2].all[pr].epitopes[e];
								lastStart[po] = deserializedDatas[2].all[pr].epitopes[e].parameters.Epitope_Start;
								/*
                                if (print[po]) // empty buffer of combined sequences.
                                {
                                    /*
                                    epitope[po].parameters = getInnovagenSolubility(epitope[po].epitopeSeq);
                                    epitope[po].parameters.Epitope_Start = lastStart[po];
                                    epitope[po].parameters.Epitope_Stop = lastEpitope[po,0].parameters.Epitope_Stop;
                                    if (!epitope[po].parameters.Innovagen_Solubility)
                                        epitope[po].parameters.Modified_Sequence_Parameters = optimize(epitope[po].epitopeSeq, "KKKK", 99, 94);
                                    if (epitope[po].parameters.Modified_Sequence_Parameters != null) //take out 5-13.
                                    {
                                        //outputs[i].Add("\t\tStart\tStop\tMolecular Weight\tExtinction Coefficient\tIso-electric point\tNet Charge at pH7\tSolubility\tGRAVY Score\tModified Sequence(If Necessary)(GRAVY Score");
                                        condense[po].Add(String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}(GRAVY: {12}, {13} Solubility)",
                                            epitope[po].epitopeSeq, epitope[po].average.ToString(),
                                            epitope[po].Binding_Core, epitope[po].parameters.Epitope_Start, epitope[po].parameters.Epitope_Stop, epitope[po].parameters.Molecular_Weight, epitope[po].parameters.Extinction_Coefficient, epitope[po].parameters.Isoelectric_point_pH,
                                            epitope[po].parameters.Net_Charge_At_pH_7, epitope[po].parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), epitope[po].parameters.GRAVY_Score, epitope[po].parameters.Modified_Sequence_Parameters.Modified_Sequence, epitope[po].parameters.Modified_Sequence_Parameters.GRAVY_Score, epitope[po].parameters.Modified_Sequence_Parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor")
                                            ));
                                    }
                                    else
                                    {
                                        condense[po].Add(String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}",
                                            epitope[po].epitopeSeq, epitope[po].average.ToString(),
                                            epitope[po].Binding_Core, 
                                            epitope[po].parameters.Epitope_Start, epitope[po].parameters.Epitope_Stop, epitope[po].parameters.Molecular_Weight, epitope[po].parameters.Extinction_Coefficient, epitope[po].parameters.Isoelectric_point_pH,
                                            epitope[po].parameters.Net_Charge_At_pH_7, epitope[po].parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), epitope[po].parameters.GRAVY_Score
                                            ));
                                    }


                                    if (ep[e].parameters.Modified_Sequence_Parameters != null)
                                    {
                                        //outputs[i].Add("\t\tStart\tStop\tMolecular Weight\tExtinction Coefficient\tIso-electric point\tNet Charge at pH7\tSolubility\tGRAVY Score\tModified Sequence(If Necessary)(GRAVY Score");
                                        condense[po].Add(String.Format("{0}\t{1}\t({2},{3},{4})\t({5},{6},{7})\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}(GRAVY: {17}, {18} Solubility)",
                                            ep[e].epitopeSeq, ep[e].populationCoverage[po].Ranking_Score.ToString(),
                                            deserializedDatas[0].all[pr].epitopes[e].Binding_Core, deserializedDatas[1].all[pr].epitopes[e].Binding_Core, deserializedDatas[2].all[pr].epitopes[e].Binding_Core,
                                            pop[po].Present_In_NetMHC2.ToString().Replace("True", "Yes").Replace("False", "No"), pop[po].Present_In_NetMHC3.ToString().Replace("True", "Yes").Replace("False", "No"), pop[po].Present_In_Tepitope.ToString().Replace("True", "Yes").Replace("False", "No"),
                                            ep[e].parameters.Epitope_Start, ep[e].parameters.Epitope_Stop, ep[e].parameters.Molecular_Weight, ep[e].parameters.Extinction_Coefficient, ep[e].parameters.Isoelectric_point_pH,
                                            ep[e].parameters.Net_Charge_At_pH_7, ep[e].parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), ep[e].parameters.GRAVY_Score, ep[e].parameters.Modified_Sequence_Parameters.Modified_Sequence, ep[e].parameters.Modified_Sequence_Parameters.GRAVY_Score, ep[e].parameters.Modified_Sequence_Parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor")
                                            ));
                                    }
                                    else
                                    {
                                        condense[po].Add(String.Format("{0}\t{1}\t({2},{3},{4})\t({5},{6},{7})\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}",
                                            ep[e].epitopeSeq, ep[e].populationCoverage[po].Ranking_Score.ToString().ToString(),
                                            deserializedDatas[0].all[pr].epitopes[e].Binding_Core, deserializedDatas[1].all[pr].epitopes[e].Binding_Core, deserializedDatas[2].all[pr].epitopes[e].Binding_Core,
                                            pop[po].Present_In_NetMHC2.ToString().Replace("True", "Yes").Replace("False", "No"), pop[po].Present_In_NetMHC3.ToString().Replace("True", "Yes").Replace("False", "No"), pop[po].Present_In_Tepitope.ToString().Replace("True", "Yes").Replace("False", "No"),
                                            ep[e].parameters.Epitope_Start, ep[e].parameters.Epitope_Stop, ep[e].parameters.Molecular_Weight, ep[e].parameters.Extinction_Coefficient, ep[e].parameters.Isoelectric_point_pH,
                                            ep[e].parameters.Net_Charge_At_pH_7, ep[e].parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), ep[e].parameters.GRAVY_Score
                                            ));
                                    }

                                    epitope[po] = new Epitope();
                                    epitope[po].epitopeSeq = deserializedDatas[2].all[pr].epitopes[e].epitopeSeq;
                                    epitope[po].average = deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Ranking_Score;
                                    epitope[po].parameters = new EpitopeParameters();
                                    epitope[po].parameters.Epitope_Stop = deserializedDatas[2].all[pr].epitopes[e].parameters.Epitope_Stop;
                                    print[po] = false;
                                    lastEpitope[po, 0] = deserializedDatas[0].all[pr].epitopes[e];
                                    lastEpitope[po, 1] = deserializedDatas[1].all[pr].epitopes[e];
                                    lastEpitope[po, 2] = deserializedDatas[2].all[pr].epitopes[e];
                                    lastStart[po] = deserializedDatas[2].all[pr].epitopes[e].parameters.Epitope_Start;
                                }
                                else
                                {
                                    /*
                                    lastStart[po] = ep[e].parameters.Epitope_Start;
                                    lastcore0[po] = deserializedDatas[0].all[pr].epitopes[e].Binding_Core;
                                    lastcore1[po] = deserializedDatas[1].all[pr].epitopes[e].Binding_Core;
                                    lastcore2[po] = deserializedDatas[2].all[pr].epitopes[e].Binding_Core;
                                    lastcoverage0[po] = deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC2;
                                    lastcoverage1[po] = deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC3;
                                    lastcoverage2[po] = deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Present_In_Tepitope;
                                    lastPosition[po] = po;
                                    combinedEpitope[po] = deserializedDatas[0].all[pr].epitopes[e].epitopeSeq;
                                    epitope[po].epitopeSeq = deserializedDatas[0].all[pr].epitopes[e].epitopeSeq;
                                    epitope[po].average = ((deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].coverage / deserializedDatas[1].all[pr].epitopes[e].populationCoverage[po].maxCoverage) / 100);
                                    
                                    epitope[po] = new Epitope();
                                    epitope[po].epitopeSeq = deserializedDatas[2].all[pr].epitopes[e].epitopeSeq;
                                    epitope[po].average = deserializedDatas[2].all[pr].epitopes[e].populationCoverage[po].Ranking_Score;
                                    epitope[po].parameters = new EpitopeParameters();
                                    epitope[po].parameters.Epitope_Stop = deserializedDatas[2].all[pr].epitopes[e].parameters.Epitope_Stop;
                                    print[po] = false;
                                    lastEpitope[po, 0] = deserializedDatas[0].all[pr].epitopes[e];
                                    lastEpitope[po, 1] = deserializedDatas[1].all[pr].epitopes[e];
                                    lastEpitope[po, 2] = deserializedDatas[2].all[pr].epitopes[e];
                                    lastStart[po] = deserializedDatas[2].all[pr].epitopes[e].parameters.Epitope_Start;

                                    if (ep[e].parameters.Modified_Sequence_Parameters != null)
                                    {
                                        //outputs[i].Add("\t\tStart\tStop\tMolecular Weight\tExtinction Coefficient\tIso-electric point\tNet Charge at pH7\tSolubility\tGRAVY Score\tModified Sequence(If Necessary)(GRAVY Score");
                                        condense[po].Add(String.Format("{0}\t{1}\t({2},{3},{4})\t({5},{6},{7})\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}(GRAVY: {17}, {18} Solubility)",
                                            ep[e].epitopeSeq, ep[e].populationCoverage[po].Ranking_Score.ToString(),
                                            deserializedDatas[0].all[pr].epitopes[e].Binding_Core, deserializedDatas[1].all[pr].epitopes[e].Binding_Core, deserializedDatas[2].all[pr].epitopes[e].Binding_Core,
                                            pop[po].Present_In_NetMHC2.ToString().Replace("True", "Yes").Replace("False", "No"), pop[po].Present_In_NetMHC3.ToString().Replace("True", "Yes").Replace("False", "No"), pop[po].Present_In_Tepitope.ToString().Replace("True", "Yes").Replace("False", "No"),
                                            ep[e].parameters.Epitope_Start, ep[e].parameters.Epitope_Stop, ep[e].parameters.Molecular_Weight, ep[e].parameters.Extinction_Coefficient, ep[e].parameters.Isoelectric_point_pH,
                                            ep[e].parameters.Net_Charge_At_pH_7, ep[e].parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), ep[e].parameters.GRAVY_Score, ep[e].parameters.Modified_Sequence_Parameters.Modified_Sequence, ep[e].parameters.Modified_Sequence_Parameters.GRAVY_Score, ep[e].parameters.Modified_Sequence_Parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor")
                                            ));
                                    }
                                    else
                                    {
                                        condense[po].Add(String.Format("{0}\t{1}\t({2},{3},{4})\t({5},{6},{7})\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}",
                                            ep[e].epitopeSeq, ep[e].populationCoverage[po].Ranking_Score.ToString().ToString(),
                                            deserializedDatas[0].all[pr].epitopes[e].Binding_Core, deserializedDatas[1].all[pr].epitopes[e].Binding_Core, deserializedDatas[2].all[pr].epitopes[e].Binding_Core,
                                            pop[po].Present_In_NetMHC2.ToString().Replace("True", "Yes").Replace("False", "No"), pop[po].Present_In_NetMHC3.ToString().Replace("True", "Yes").Replace("False", "No"), pop[po].Present_In_Tepitope.ToString().Replace("True", "Yes").Replace("False", "No"),
                                            ep[e].parameters.Epitope_Start, ep[e].parameters.Epitope_Stop, ep[e].parameters.Molecular_Weight, ep[e].parameters.Extinction_Coefficient, ep[e].parameters.Isoelectric_point_pH,
                                            ep[e].parameters.Net_Charge_At_pH_7, ep[e].parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), ep[e].parameters.GRAVY_Score
                                            ));
                                    }
                                }*/
							}


							//outputs[po].Add(deserializedDatas[0].all[pr].epitopes[e].epitopeSeq +"\t" + deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC2.ToString() + "\t" + deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].Present_In_NetMHC3.ToString() + "\t" + deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].Present_In_Tepitope.ToString() + "\t" + deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].coverage.ToString() + "\t" + deserializedDatas[0].all[pr].epitopes[e].populationCoverage[po].maxCoverage.ToString());
						}
					}
				}

				for (int i = 0; i < outputs.Length; i++)
				{
					//condense[i].Add("");
					outputs[i].Add("");
				}
			}

			for (int i = 0; i < outputs.Length; i++)
			{
				l.write(outputs[i], path: resultDir + "/Result_Summary/Raw_" + outputs[i][0].Replace("#", ""), delete: false);
				//l.write(condense[i], path: resultDir + "/Result_Summary/Summarized_" + condense[i][0].Replace("#", ""), delete: false);
			}

			foreach (string s in Directory.GetFiles(resultDir + "/Result_Summary"))
			{
				if (s.Contains("Raw_") && !s.Contains(".f."))
				{
					//Epitope Sequence	Overall Score	Binding Core(NetMHC2, NetMHC3, Tepitope)	Predicted Score(NetMHC2, NetMHC3, Tepitope)	Population Coverage(NetMHC2/Max Coverage, NetMHC3/Max Coverage, Tepitope/Max Coverage)	Predicted by(NetMHC2, NetMHC3, Tepitope)	Start	Stop	Molecular Weight	Extinction Coefficient	Iso-electric point	Net Charge at pH7	Solubility	GRAVY Score	Modified Sequence	Modified GRAVY Score	Modified Solubility)
					List<string> input = l.read(path: s);
					List<string> output2 = new List<string>();
					int lastStart2 = -99;
					double lastScore2 = 999999;
					string core = "";
					foreach (string i in input)
					{
						if (i.Contains("Antigen") || i.Contains("Epitope") || i.Contains("#") || i == "")
						{
							output2.Add(i);
						}
						else
						{
							if (int.Parse(returnPositionText(i, 6)) == lastStart2 + 1)
							{
								if (double.Parse(returnPositionText(i, 1)) < lastScore2) //if score less than the last
								{
									output2.RemoveAt(output2.Count - 1);
									output2.Add(i); //.Replace(returnPositionText(i, 2), core));
									lastScore2 = double.Parse(returnPositionText(i, 1));
									//if (!core.Contains(returnPositionText(i, 2).Substring(1).Substring(0, returnPositionText(i, 2).Substring(1).IndexOf(","))))
									// core += returnPositionText(i, 2).Substring(1).Substring(0, returnPositionText(i, 2).Substring(1).IndexOf(",")) + ",";
								}
								lastStart2 = int.Parse(returnPositionText(i, 6));
							}
							else
							{
								//core = returnPositionText(i, 2).Substring(1).Substring(0, returnPositionText(i, 2).Substring(1).IndexOf(",")) + ",";
								output2.Add(i);//.Replace(returnPositionText(i, 2), core));
								lastStart2 = int.Parse(returnPositionText(i, 6));
								lastScore2 = double.Parse(returnPositionText(i, 1));
							}
						}
					}

					for (int str = 0; str < output2.Count; str++)
						output2[str] = output2[str].Replace(",\t", "\t");
					l.write(output2, path: s + ".f.txt");
				}
			}
		}
		//rank the summarized results and pull some parameter data
		public void rankResult()
		{
			foreach (string s in Directory.GetFiles(resultDir + "/Result_Summary"))
			{
				if (s.Contains(".f.")) //s.Contains("Summarized_") || 
				{
					List<string> output = new List<string>();
					Dictionary<double, string> list = new Dictionary<double, string>();
					EpitopeParameters ep = new EpitopeParameters();
					int counter = 1;
					foreach (string t in l.read(path: s))
					{
						if (t.Contains("#"))
						{
							output.Add(t);
						}
						else if (t.Contains("Antigen"))
						{
							output.Add(t);
						}
						else if (t.Contains("Epitope Sequence"))
						{
							output.Add("Rank\t" + t);
							//labels
						}
						else if (t == "")
						{
							//break between antigen
							list = list.OrderBy(x => x.Key).ToDictionary(x => x.Key, x => x.Value);
							foreach (KeyValuePair<double, string> k in list)
							{
								output.Add(counter.ToString() + "\t" + k.Value);
								counter++;
							}
							list = new Dictionary<double, string>();
							counter = 1;
							ep = new EpitopeParameters();
						}
						else
						{
							//do something with the binding cores
							if (t.Split('\t').Length - 1 < 7)//!@#!@#!@#!@#! CHANGE TO 7
							{
								//need to find parameter
								ep = getInnovagenSolubility(t.Substring(0, t.IndexOf("\t")));
								//Molecular Weight	Extinction Coefficient	Iso-electric point	Net Charge at pH7	Solubility	GRAVY Score	Modified Sequence(GRAVY Score, Solubility)
								if (ep.Innovagen_Solubility)
								{
									try
									{
										list.Add(double.Parse(returnPositionText(t, 1)), string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}", t, ep.Molecular_Weight, ep.Extinction_Coefficient, ep.Isoelectric_point_pH, ep.Net_Charge_At_pH_7, ep.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), ep.GRAVY_Score));
									}
									catch
									{

									}
								}
								else
								{
									try
									{
										ep.Modified_Sequence_Parameters = optimize(t.Substring(0, t.IndexOf("\t")), "KKKK", 99, 94);
										list.Add(double.Parse(returnPositionText(t, 1)), string.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9} Solubility)", t, ep.Molecular_Weight, ep.Extinction_Coefficient, ep.Isoelectric_point_pH, ep.Net_Charge_At_pH_7, ep.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor"), ep.GRAVY_Score, ep.Modified_Sequence_Parameters.Modified_Sequence, ep.Modified_Sequence_Parameters.GRAVY_Score, ep.Modified_Sequence_Parameters.Innovagen_Solubility.ToString().Replace("True", "Good").Replace("False", "Poor")));
									}
									catch
									{

									}
								}
							}
							else
								try
							{
								list.Add(double.Parse(returnPositionText(t, 1)), t);
							}

							catch
							{

							}
						}
					}

					l.write(output, path: s, delete: true);
				}
			}
		}
		//QoL function used to return a substring in a tab delimitted file
		public string returnPositionText(string s, int tabNum)
		{
			string i = "";
			try
			{
				i = s.Substring(l.NthOccurence(s, Convert.ToChar("\t"), tabNum) + 1, l.NthOccurence(s, Convert.ToChar("\t"), tabNum + 1) - l.NthOccurence(s, Convert.ToChar("\t"), tabNum) - 1);
			}
			catch
			{
				throw new Exception("0j908");
			}
			return i;
		}
		#endregion

		/// <summary>
		/// function used to generation the final html report
		/// </summary>
		#region generate the "readable" report
		int figure_index = 1;
		int table_index = 1;
		int epitope_index = 1;
		public void generateReport()
		{
			List<string> output = new List<string>();
			//List<string> template = l.read(path: resultDir + "/template.html");
			List<Protein> proteins = new List<Protein>();
			Protein p = new Protein();
			List<Epitope> epitopes = new List<Epitope>();

			foreach (string f in Directory.GetFiles(resultDir + "/Result_Summary"))
			{
				if (f.Contains("Summarized_world"))
				{
					List<string> input = l.read(path: f);
					bool b = false;
					for(int str = 0; str < input.Count; str ++)
					{
						if (input[str].Contains("#"))
						{

						}
						else if (input[str].Contains("Epitope Sequence"))
						{

						}
						else if (input[str].Contains("Antigen") || input[str].Contains(".fasta"))
						{
							if (b)
							{
								p.epitopes = epitopes;
								proteins.Add(p);
								p = new Protein();
								epitopes = new List<Epitope>();
							}
							//proteins.proteinName = input[str].Substring(input[str].IndexOf("."));
							p.proteinName = input[str].Substring(input[str].IndexOf(": ") + 2);
							b = true;

						}
						else if (input[str] == "" || str == input.Count - 1)
						{
							p.epitopes = epitopes;
							proteins.Add(p);
							p = new Protein();
							epitopes = new List<Epitope>();
						}
						else
						{
							Epitope e = new Epitope();
							//Epitope Sequence	Overall Score	Binding Core	Start	Stop	Molecular Weight	Extinction Coefficient	Iso-electric point	Net Charge at pH7	Solubility	GRAVY Score	Modified Sequence(GRAVY Score, Solubility)
							e.epitopeSeq = returnPositionText(input[str], 1);
							e.average = double.Parse(returnPositionText(input[str], 2));
							e.Binding_Core = returnPositionText(input[str], 3).Split(Convert.ToChar(",")).ToList();
							e.parameters = new EpitopeParameters();
							e.parameters.Epitope_Start = int.Parse(returnPositionText(input[str], 4));
							e.parameters.Isoelectric_point_pH = double.Parse(returnPositionText(input[str], 8));
							e.parameters.GRAVY_Score = double.Parse(returnPositionText(input[str], 11));
							try {
								e.parameters.Modified_Sequence_Parameters = new ModifiedParameters();
								e.parameters.Modified_Sequence_Parameters.Modified_Sequence = returnPositionText(input[str], 12);
								e.parameters.Modified_Sequence_Parameters.GRAVY_Score = double.Parse(returnPositionText(input[str], 13));
							}
							catch
							{
								e.parameters.Modified_Sequence_Parameters = null;
							}
							epitopes.Add(e);
						}
					}

				}
			}
			output.Add("<!DOCTYPE html> <html lang=\"en\"> <style> table, th, td {border: 1px solid black; border-collapse: collapse; } th, td { padding: 5px; } th { text-align: center; } </style> <body>"); //add header
			foreach (Protein pr in proteins)
			{
				List<string> templateMod = new List<string>();
				templateMod = l.read(path: resultDir + "/template.html");
				int i = 0;
				int max = templateMod.Count;
				while (i < max)
				{
					while (templateMod[i].Contains("{") && templateMod[i].Contains("}"))
					{
						List<string> input;
						string s = templateMod[i].Substring(templateMod[i].IndexOf("{") + 1, templateMod[i].IndexOf("}") - templateMod[i].IndexOf("{") - 1);
						switch (s)
						{
							case "antigen":
							templateMod[i] = templateMod[i].Replace("{antigen}", pr.proteinName);
							break;
							case "index":
							templateMod[i] = templateMod[i].Replace("{index}", epitope_index.ToString());
							epitope_index++;
							break;
							case "antigen_info":
							templateMod.RemoveAt(i);
							input = l.read(path: root + "/antigen_info.txt");
							templateMod.InsertRange(i, input);
							i += input.Count - 1;
							max = templateMod.Count;
							break;
							case "fig_index":
							templateMod[i] = templateMod[i].Replace("{fig_index}", figure_index.ToString());
							figure_index++;
							break;
							case "img_id":
							templateMod[i] = templateMod[i].Replace("{img_id}", pr.proteinName);
							break;
							case "data_fields":
							templateMod.RemoveAt(i);
							for (int ep = 0; ep < pr.epitopes.Count; ep++)// ep in pr.epitopes)
							{
								List<string> temp = new List<string>();
								temp.Add("<tr>");
								temp.Add("<th>" + pr.proteinName + "_" + (ep+1).ToString() + "</th>");
								temp.Add("<th>" + pr.epitopes[ep].average.ToString() + "</th>");
								temp.Add("<th>" + pr.epitopes[ep].epitopeSeq + "</th>");
								foreach (string bc in pr.epitopes[ep].Binding_Core)
									temp.Add("<th>" + bc.Replace(",","<br>") + "</th>");
								temp.Add("<th>" + pr.epitopes[ep].parameters.Epitope_Start.ToString() + "</th>");
								temp.Add("<th>" + pr.epitopes[ep].parameters.GRAVY_Score.ToString() + "</th>");
								temp.Add("<th>" + pr.epitopes[ep].parameters.Isoelectric_point_pH.ToString() + "</th>");
								if (pr.epitopes[ep].parameters.Modified_Sequence_Parameters == null)
									temp.Add("<th>" + "Soluble" + "</th>");
								else
									temp.Add("<th>" + pr.epitopes[ep].parameters.Modified_Sequence_Parameters.Modified_Sequence + "(" + pr.epitopes[ep].parameters.Modified_Sequence_Parameters.GRAVY_Score + ") </th>");
								temp.Add("</tr>");
								templateMod.InsertRange(i, temp);
								i += temp.Count - 1;
								max = templateMod.Count;
							}

							break;
							case "table_index":
							templateMod[i] = templateMod[i].Replace("{table_index}", table_index.ToString());
							table_index++;
							break;
							case "prev_pub_description":
							templateMod.RemoveAt(i);
							input = l.read(path: root + "/prev_epitopes.txt");
							templateMod.InsertRange(i, input);
							i += input.Count - 1;
							max = templateMod.Count;
							break;
							case "common_fields":
							templateMod[i] = templateMod[i].Replace("{common_fields}", "");
							/*List<string> temp1 = new List<string>();
                                foreach (string str in l.read(root+"/common.txt"))
                                {
                                    bool start = false;
                                    if (start)
                                    {
                                        if (!str.Contains("#") && str != "")
                                        {
                                            temp1.Add(str);
                                        }
                                        else
                                            start = false;
                                    }
                                    if (s.Contains("#" + pr.proteinName))
                                        start = true;
                                }
                                foreach (string knownEpitope in temp1)
                                {

                                }*/
							break;
							case "conclusion":
							templateMod.RemoveAt(i);
							input = l.read(path: root + "/conclusion.txt");
							templateMod.InsertRange(i, input);
							i += input.Count - 1;
							max = templateMod.Count;
							break;
							case "table_sub_index":
							break;
							case "Populations":
							break;
							default:
							throw new Exception("0j606 - Unknown placeholder used: " + s);
						}
					}
					i++;
				}
				output = output.Concat(templateMod).ToList();
				//templateMod = template;
			}

			output.Add("</body>");
			l.write(output, path: resultDir + "/test.html");
		}
		#endregion

		/// <summary>
		/// function used to modify the peptide to make it soluble
		/// </summary>
		/// <param name="epitope"></param>
		/// <param name="protein">whole antigen sequence</param>
		/// <param name="start"></param>
		/// <param name="stop"></param>
		/// <param name="innovagen"></param>
		/// <returns></returns>
		#region calculate the epitope parameters
		public ModifiedParameters optimize(string peptide, string protein, int start, int stop)
		{
			string epitope = peptide.Replace("C", "S");
			string modifiedPeptide = "!err";
			string tempPep = "!err";
			string longest = "!err";
			try
			{
				longest = protein.Substring(start - 3, (stop - start) + 6);
			}
			catch (Exception e)
			{
				longest = "KKK" + epitope + "KKK";
			}
			bool soluble = false;

			char[] left = longest.Substring(0, 3).ToCharArray();
			char[] right = longest.Substring(longest.Length - 3, 3).ToCharArray();

			for (int leftCount = left.Length - 1; leftCount >= 0; leftCount--)
			{ //add 1 to left side
				if (leftCount == left.Length - 1)
				{
					modifiedPeptide = left[leftCount] + epitope;
					tempPep = modifiedPeptide;
				}
				else
				{
					modifiedPeptide = left[leftCount] + tempPep;
					tempPep = modifiedPeptide;
				}
				soluble = getInnovagenSolubility(modifiedPeptide).Innovagen_Solubility;

				if (soluble) break;

				for (int rightCount = 0; rightCount < right.Length; rightCount++)
				{ //add 1 to right side
					modifiedPeptide += right[rightCount];
					soluble = getInnovagenSolubility(modifiedPeptide).Innovagen_Solubility;

					if (soluble) break;
					//6
				}
				if (soluble) break;
			}

			if (!soluble)
			{
				for (int leftCount = 3; leftCount > 0; leftCount--)
				{ //add 1 to left side
					modifiedPeptide = "K" + epitope;

					soluble = getInnovagenSolubility(modifiedPeptide).Innovagen_Solubility;
					if (soluble) break;

					for (int rightCount = 0; rightCount < 3; rightCount++)
					{ //add 1 to right side
						modifiedPeptide += "K";
						soluble = getInnovagenSolubility(modifiedPeptide).Innovagen_Solubility;
						if (soluble) break;
						//6
					}
					if (soluble) break;
				}
			}

			if (modifiedPeptide == "!err")
				throw new Exception(string.Format("error 0j107 - error modifying peptide: {0}, {1}, {2}, {3}, {4}, {5}", epitope, protein, start, stop, modifiedPeptide));

			ModifiedParameters m = new ModifiedParameters();
			if (!soluble)
			{
				m.Modified_Sequence = "Unable to find soluble modification";
				m.Innovagen_Solubility = false;
				m.GRAVY_Score = -9999999;
			}
			else
			{
				EpitopeParameters t = getInnovagenSolubility(modifiedPeptide);
				m.Modified_Sequence = modifiedPeptide;
				m.Isoelectric_point_pH = t.Isoelectric_point_pH;
				m.Net_Charge_At_pH_7 = t.Net_Charge_At_pH_7;
				m.Innovagen_Solubility = t.Innovagen_Solubility;
				m.GRAVY_Score = calculateGRAVY(modifiedPeptide);
				if (Directory.Exists (root + "/Modified_Protein_Sequence")) 
				{
					Directory.Delete (root + "/Modified_Protein_Sequence",true);
					Directory.CreateDirectory (root + "/Modified_Protein_Sequence");
				}
				System.IO.File.WriteAllText(root+"/Modified_Protein_Sequence/" + m.Modified_Sequence, m.Modified_Sequence);
				callEntryPeptide ();
				List<parsemodified> p = parseModifiedResults ();
				if (p.Count > 1)
					throw new Exception ("more than 1 modified sequence??? how");
				m.Alleles = new List<Allele> ();
				List<double> sum = new List<double> ();
				foreach (KeyValuePair<string, double> k in p[0].alleleAffinity) 
				{
					Allele a = new Allele ();
					a.alleleName = k.Key;
					a.IC50 = k.Value;
					m.Alleles.Add (a);
					sum.Add(a.IC50);
					a = new Allele ();
				}
				m.Modified_Average = trimmedMean (sum, 0.2);
				sum = new List<double> ();
				m.Modified_Cores = new List<string> ();
				foreach (KeyValuePair<string,string> k in p[0].alleleCore) 
				{
					if (!m.Modified_Cores.Contains(k.Value))
						m.Modified_Cores.Add (k.Value);
				}
			}

			return m;
		}
		//get innovagen solubility, requires dll ver 1231 or higher.
		public EpitopeParameters getInnovagenSolubility(string peptide)
		{
			//parameters used to query innovagen predictor.
			string URI = "http://pepcalc.com/ppc.php";
			Dictionary<string, string> reqparm = new Dictionary<string, string>();
			reqparm.Add("sequence", peptide);
			reqparm.Add("nTerm", "(NH2-)");
			reqparm.Add("cTerm", "(-COOH)");
			reqparm.Add("aaCode", "0");

			//get result
			List<string> phpResult = l.getDataFromPHP(URI, "POST", reqparm);

			//parse the php data stream to get necessary data.
			string weight = "!err";
			string extinction = "!err";
			string pi = "!err";
			string charge = "!err";
			string solubility = "!err";

			for (int i = 0; i < phpResult.Count; i++)
			{
				if (phpResult[i].Contains("Molecular weight:"))
					weight = phpResult[i + 1].Substring(phpResult[i + 1].IndexOf(">") + 1, phpResult[i + 1].IndexOf("</td>") - phpResult[i + 1].IndexOf(">") - 1);
				else if (phpResult[i].Contains("Extinction coefficient:"))
					extinction = phpResult[i + 1].Substring(phpResult[i + 1].IndexOf(">") + 1, phpResult[i + 1].IndexOf("</td>") - phpResult[i + 1].IndexOf(">") - 1);
				else if (phpResult[i].Contains("Iso-electric point:"))
					pi = phpResult[i + 1].Substring(phpResult[i + 1].IndexOf(">") + 1, phpResult[i + 1].IndexOf("</td>") - phpResult[i + 1].IndexOf(">") - 1);
				else if (phpResult[i].Contains("Net charge at pH 7:"))
					charge = phpResult[i + 1].Substring(phpResult[i + 1].IndexOf(">") + 1, phpResult[i + 1].IndexOf("</td>") - phpResult[i + 1].IndexOf(">") - 1);
				else if (phpResult[i].Contains("Estimated solubility:"))
					solubility = phpResult[i + 1].Substring(phpResult[i + 1].IndexOf(">") + 1, 4);
			}

			if (weight == "!err" || extinction == "!err" || pi == "!err" || charge == "!err" || solubility == "!err")
				throw new Exception(string.Format("error 0j101 - unexpected character in getInnovagenSolubility {0}, {1}, {2}, {3}, {4}", weight, extinction, pi, charge, solubility));

			EpitopeParameters e = new EpitopeParameters();
			string temp = weight.Substring(0, weight.IndexOf(" g"));
			e.Molecular_Weight = double.Parse(temp);
			temp = extinction.Substring(0, extinction.IndexOf(" M"));
			e.Extinction_Coefficient = int.Parse(temp);
			temp = pi.Substring(3);
			e.Isoelectric_point_pH = double.Parse(temp);
			e.Net_Charge_At_pH_7 = double.Parse(charge);
			e.GRAVY_Score = calculateGRAVY(peptide);
			if (solubility == "Good")
				e.Innovagen_Solubility = true;
			else
				e.Innovagen_Solubility = false;
			return e;
		}
		public double calculateGRAVY(string _seq)
		{
			Dictionary<string, double> hashes = getGravyHash();
			string seq = _seq.ToLower();

			double gravyResult = 0;
			//The GRAVY value for a peptide or protein is calculated as the sum of hydropathy values [9]
			//of all the amino acids, divided by the number of residues in the sequence. 

			for (var i = 0; i < seq.Length; i++)
			{
				gravyResult = gravyResult + hashes[seq[i].ToString()];
			}
			if (seq.Length > 0)
			{
				gravyResult = gravyResult / seq.Length;
			}
			else
			{
				Console.WriteLine("The sequence is too short");
			}
			return Math.Round(gravyResult, 2);
		}
		private Dictionary<string, double> getGravyHash()
		{
			//Author(s): Kyte J., Doolittle R.F.
			//Reference: J. Mol. Biol. 157:105-132(1982).	
			Dictionary<string, double> hash = new Dictionary<string, double>();
			hash.Add("a", 1.800);
			hash.Add("r", -4.500);
			hash.Add("n", -3.500);
			hash.Add("d", -3.500);
			hash.Add("c", 2.500);
			hash.Add("q", -3.500);
			hash.Add("e", -3.500);
			hash.Add("g", -0.400);
			hash.Add("h", -3.200);
			hash.Add("i", 4.500);
			hash.Add("l", 3.800);
			hash.Add("k", -3.900);
			hash.Add("m", 1.900);
			hash.Add("f", 2.800);
			hash.Add("p", -1.600);
			hash.Add("s", -0.800);
			hash.Add("t", -0.700);
			hash.Add("w", -0.900);
			hash.Add("y", -1.300);
			hash.Add("v", 4.200);
			return hash;
		}
		#endregion

		#region writting debug info
		private void writeDebug()
		{
			File.Delete(root + "/Debug_Settings.json");
			using (StreamWriter f = File.AppendText(root + "/Debug_Settings.json"))
			{
				JsonSerializer serializer = new JsonSerializer();
				serializer.Serialize(f, DS);
			}
		}
		private void readDebug()
		{
			using (StreamReader f = File.OpenText(root + "/Debug_Settings.json"))
			{
				JsonSerializer serializer = new JsonSerializer();
				DS = (DebugSettings)serializer.Deserialize(f, typeof(DebugSettings));
			}
		}
		class DebugSettings
		{
			public int stepsCompleted { get; set; }
			public double cutoff { get; set; }
			public double analysisTime { get; set; }
			public bool isError { get; set; }
			public bool nonPubSequence { get; set; }
		}
		#endregion

		//stuff to make figure
		public void makeBexTable()
		{
			//table format: Name > Allegen > Epitope > Cores (2.2) > start > Tepitope score.
			//input: world.f.txt, all data is in there. 
			rawData bexData = new rawData();
			rawData bexData2 = new rawData();
			Protein bexProtein = new Protein();
			Epitope bexEpitope = new Epitope();
			List<Protein> bexProtList = new List<Protein>();
			List<Protein> bexProtList2 = new List<Protein>();
			bexData.all = bexProtList;
			bexData2.all = bexProtList2;
			List<Epitope> bexEpList = new List<Epitope>();
			Allele bexAllele = new Allele();
			List<Allele> bexAlleleList = new List<Allele>();
			bool start = false;
			List<string> input = l.read(path: resultDir + "/Result_Summary" + "/Raw_world.f.txt");

			#region readjson
			readJsonIntoMemory();
			/*
            foreach (string s in Directory.GetFiles(resultDir)) //read in the 3 json files.
            {
                if (s.Contains(".json"))
                {
                    if (s.Contains("netMHCII-2.2_results.json"))
                    {
                        using (StreamReader f = File.OpenText(s))
                        {
                            JsonSerializer serializer = new JsonSerializer();
                            deserializedDatas[0] = (rawData)serializer.Deserialize(f, typeof(rawData));
                        }
                    }
                    else if (s.Contains("iedb_results.json"))
                    {
                        using (StreamReader f = File.OpenText(s))
                        {
                            JsonSerializer serializer = new JsonSerializer();
                            deserializedDatas[2] = (rawData)serializer.Deserialize(f, typeof(rawData));
                        }
                    }
                }
            }*/
			#endregion
			#region loadData
			//0 = 2.2 = bexdata
			foreach (string s in input)
			{
				if (s.Contains("#"))
				{
					//ignore
				}
				else if (s.Contains("Antigen:"))
				{
					if (!start) start = true;
					else
					{
						bexProtein.epitopes = bexEpList;
						bexEpList = new List<Epitope>();
						bexData.all.Add(bexProtein);
					}
					bexProtein = new Protein();
					bexProtein.proteinName = s.Substring(s.IndexOf(": ") + 2);
				}
				else if (s.Contains("Rank"))
				{
					//igore too
				}
				else if (input.LastIndexOf(s) == input.Count - 1)
				{
					//   Rank Epitope Sequence Overall Score Binding Core(NetMHC2, NetMHC3, Tepitope)    Predicted Score(NetMHC2, NetMHC3, Tepitope) Population Coverage(NetMHC2 / Max Coverage, NetMHC3 / Max Coverage, Tepitope / Max Coverage)  Predicted by(NetMHC2, NetMHC3, Tepitope)    Start   Stop    Molecular Weight    Extinction Coefficient  Iso - electric point  Net Charge at pH7   Solubility  GRAVY Score Modified Sequence   Modified GRAVY Score    Modified Solubility)
					//  1   MGISFILAHTYGYPR 25.67   YRMGISFIL,ISFILAHTY(197.99, 142.31, 2.80)(52.91002 / 57.49, 88.05464 / 100.02, 60.66020 / 69.40)(Yes, No, No) 324 338 1726.01 2560    9.39    1.1 Poor    0.34    KMGISFILAHTYGYPRKK - 0.37   Good
					bexEpitope.epitopeSeq = returnPositionText(s, 1);
					bexEpitope.average = double.Parse(returnPositionText(s, 2));
					bexEpitope.Binding_Core = returnPositionText(s, 3).Substring(1, returnPositionText(s, 3).IndexOf("|") - 1).Split(Convert.ToChar(",")).ToList();
					bexEpitope.parameters = new EpitopeParameters();
					bexEpitope.parameters.Epitope_Start = int.Parse(returnPositionText(s, 7));
					bexEpitope.parameters.Epitope_Stop = int.Parse(returnPositionText(s, 8));
					bexEpitope.parameters.Molecular_Weight = double.Parse(returnPositionText(s, 9));
					bexEpitope.parameters.Isoelectric_point_pH = double.Parse(returnPositionText(s, 11));
					bexEpitope.parameters.Net_Charge_At_pH_7 = double.Parse(returnPositionText(s, 12));
					if (returnPositionText(s, 13) == "Good")
					{
						bexEpitope.parameters.Innovagen_Solubility = true;
						bexEpitope.parameters.GRAVY_Score = double.Parse(s.Substring(s.LastIndexOf("\t") + 1));
					}
					else
					{
						bexEpitope.parameters.GRAVY_Score = double.Parse(returnPositionText(s, 14));
						bexEpitope.parameters.Innovagen_Solubility = true;
						bexEpitope.parameters.Modified_Sequence_Parameters = new ModifiedParameters();
						bexEpitope.parameters.Modified_Sequence_Parameters.Modified_Sequence = returnPositionText(s, 15);
						bexEpitope.parameters.Modified_Sequence_Parameters.GRAVY_Score = double.Parse(returnPositionText(s, 16));
					}
					//find the allele data
					foreach (Protein p in deserializedDatas[0].all)
					{
						if (p.proteinName == bexProtein.proteinName)
						{
							foreach (Epitope e in p.epitopes)
							{
								if (e.epitopeSeq == bexEpitope.epitopeSeq)
								{
									bexEpitope.alleles = e.alleles;
								}
							}
						}
					}

					bexEpList.Add(bexEpitope);
					bexEpitope = new Epitope();
					bexProtein.epitopes = bexEpList;
					bexEpList = new List<Epitope>();
					bexData.all.Add(bexProtein);

					bexProtein = new Protein();
				}
				else
				{
					//   Rank Epitope Sequence Overall Score Binding Core(NetMHC2, NetMHC3, Tepitope)    Predicted Score(NetMHC2, NetMHC3, Tepitope) Population Coverage(NetMHC2 / Max Coverage, NetMHC3 / Max Coverage, Tepitope / Max Coverage)  Predicted by(NetMHC2, NetMHC3, Tepitope)    Start   Stop    Molecular Weight    Extinction Coefficient  Iso - electric point  Net Charge at pH7   Solubility  GRAVY Score Modified Sequence   Modified GRAVY Score    Modified Solubility)
					//  1   MGISFILAHTYGYPR 25.67   YRMGISFIL,ISFILAHTY(197.99, 142.31, 2.80)(52.91002 / 57.49, 88.05464 / 100.02, 60.66020 / 69.40)(Yes, No, No) 324 338 1726.01 2560    9.39    1.1 Poor    0.34    KMGISFILAHTYGYPRKK - 0.37   Good
					bexEpitope.epitopeSeq = returnPositionText(s, 1);
					bexEpitope.average = double.Parse(returnPositionText(s, 2));
					bexEpitope.Binding_Core = returnPositionText(s, 3).Substring(1, returnPositionText(s, 3).IndexOf("|") - 1).Split(Convert.ToChar(",")).ToList();
					bexEpitope.parameters = new EpitopeParameters();
					bexEpitope.parameters.Epitope_Start = int.Parse(returnPositionText(s, 7));
					bexEpitope.parameters.Epitope_Stop = int.Parse(returnPositionText(s, 8));
					bexEpitope.parameters.Molecular_Weight = double.Parse(returnPositionText(s, 9));
					bexEpitope.parameters.Isoelectric_point_pH = double.Parse(returnPositionText(s, 11));
					bexEpitope.parameters.Net_Charge_At_pH_7 = double.Parse(returnPositionText(s, 12));
					if (returnPositionText(s, 13) == "Good")
					{
						bexEpitope.parameters.Innovagen_Solubility = true;
						bexEpitope.parameters.GRAVY_Score = double.Parse(s.Substring(s.LastIndexOf("\t") + 1));
					}
					else
					{
						bexEpitope.parameters.GRAVY_Score = double.Parse(returnPositionText(s, 14));
						bexEpitope.parameters.Innovagen_Solubility = true;
						bexEpitope.parameters.Modified_Sequence_Parameters = new ModifiedParameters();
						bexEpitope.parameters.Modified_Sequence_Parameters.Modified_Sequence = returnPositionText(s, 15);
						bexEpitope.parameters.Modified_Sequence_Parameters.GRAVY_Score = double.Parse(returnPositionText(s, 16));
					}
					//find the allele data
					foreach (Protein p in deserializedDatas[0].all)
					{
						if (p.proteinName == bexProtein.proteinName)
						{
							foreach (Epitope e in p.epitopes)
							{
								if (e.epitopeSeq == bexEpitope.epitopeSeq)
								{
									bexEpitope.alleles = e.alleles;
								}
							}
						}
					}

					bexEpList.Add(bexEpitope);
					bexEpitope = new Epitope();
				}
			}

			//2 = iedb = bexdata2
			start = false;
			foreach (string s in input)
			{
				if (s.Contains("#"))
				{
					//ignore
				}
				else if (s.Contains("Antigen:"))
				{
					if (!start) start = true;
					else
					{
						bexProtein.epitopes = bexEpList;
						bexEpList = new List<Epitope>();
						bexData2.all.Add(bexProtein);
					}
					bexProtein = new Protein();
					bexProtein.proteinName = s.Substring(s.IndexOf(": ") + 2);
				}
				else if (s.Contains("Rank"))
				{
					//igore too
				}
				else if (input.LastIndexOf(s) == input.Count - 1)
				{
					//   Rank Epitope Sequence Overall Score Binding Core(NetMHC2, NetMHC3, Tepitope)    Predicted Score(NetMHC2, NetMHC3, Tepitope) Population Coverage(NetMHC2 / Max Coverage, NetMHC3 / Max Coverage, Tepitope / Max Coverage)  Predicted by(NetMHC2, NetMHC3, Tepitope)    Start   Stop    Molecular Weight    Extinction Coefficient  Iso - electric point  Net Charge at pH7   Solubility  GRAVY Score Modified Sequence   Modified GRAVY Score    Modified Solubility)
					//  1   MGISFILAHTYGYPR 25.67   YRMGISFIL,ISFILAHTY(197.99, 142.31, 2.80)(52.91002 / 57.49, 88.05464 / 100.02, 60.66020 / 69.40)(Yes, No, No) 324 338 1726.01 2560    9.39    1.1 Poor    0.34    KMGISFILAHTYGYPRKK - 0.37   Good
					bexEpitope.epitopeSeq = returnPositionText(s, 1);
					bexEpitope.average = double.Parse(returnPositionText(s, 2));
					bexEpitope.Binding_Core = returnPositionText(s, 3).Substring(1, returnPositionText(s, 3).IndexOf("|") - 1).Split(Convert.ToChar(",")).ToList();
					bexEpitope.parameters = new EpitopeParameters();
					bexEpitope.parameters.Epitope_Start = int.Parse(returnPositionText(s, 7));
					bexEpitope.parameters.Epitope_Stop = int.Parse(returnPositionText(s, 8));
					bexEpitope.parameters.Molecular_Weight = double.Parse(returnPositionText(s, 9));
					bexEpitope.parameters.Isoelectric_point_pH = double.Parse(returnPositionText(s, 11));
					bexEpitope.parameters.Net_Charge_At_pH_7 = double.Parse(returnPositionText(s, 12));
					if (returnPositionText(s, 13) == "Good")
					{
						bexEpitope.parameters.Innovagen_Solubility = true;
						bexEpitope.parameters.GRAVY_Score = double.Parse(s.Substring(s.LastIndexOf("\t") + 1));
					}
					else
					{
						bexEpitope.parameters.GRAVY_Score = double.Parse(returnPositionText(s, 14));
						bexEpitope.parameters.Innovagen_Solubility = true;
						bexEpitope.parameters.Modified_Sequence_Parameters = new ModifiedParameters();
						bexEpitope.parameters.Modified_Sequence_Parameters.Modified_Sequence = returnPositionText(s, 15);
						bexEpitope.parameters.Modified_Sequence_Parameters.GRAVY_Score = double.Parse(returnPositionText(s, 16));
					}
					//find the allele data
					foreach (Protein p in deserializedDatas[2].all)
					{
						if (p.proteinName == bexProtein.proteinName)
						{
							foreach (Epitope e in p.epitopes)
							{
								if (e.epitopeSeq == bexEpitope.epitopeSeq)
								{
									bexEpitope.alleles = e.alleles;
								}
							}
						}
					}

					bexEpList.Add(bexEpitope);
					bexEpitope = new Epitope();
					bexProtein.epitopes = bexEpList;
					bexEpList = new List<Epitope>();
					bexData2.all.Add(bexProtein);

					bexProtein = new Protein();
				}
				else
				{
					//   Rank Epitope Sequence Overall Score Binding Core(NetMHC2, NetMHC3, Tepitope)    Predicted Score(NetMHC2, NetMHC3, Tepitope) Population Coverage(NetMHC2 / Max Coverage, NetMHC3 / Max Coverage, Tepitope / Max Coverage)  Predicted by(NetMHC2, NetMHC3, Tepitope)    Start   Stop    Molecular Weight    Extinction Coefficient  Iso - electric point  Net Charge at pH7   Solubility  GRAVY Score Modified Sequence   Modified GRAVY Score    Modified Solubility)
					//  1   MGISFILAHTYGYPR 25.67   YRMGISFIL,ISFILAHTY(197.99, 142.31, 2.80)(52.91002 / 57.49, 88.05464 / 100.02, 60.66020 / 69.40)(Yes, No, No) 324 338 1726.01 2560    9.39    1.1 Poor    0.34    KMGISFILAHTYGYPRKK - 0.37   Good
					bexEpitope.epitopeSeq = returnPositionText(s, 1);
					bexEpitope.average = double.Parse(returnPositionText(s, 2));
					bexEpitope.Binding_Core = returnPositionText(s, 3).Substring(1, returnPositionText(s, 3).IndexOf("|") - 1).Split(Convert.ToChar(",")).ToList();
					bexEpitope.parameters = new EpitopeParameters();
					bexEpitope.parameters.Epitope_Start = int.Parse(returnPositionText(s, 7));
					bexEpitope.parameters.Epitope_Stop = int.Parse(returnPositionText(s, 8));
					bexEpitope.parameters.Molecular_Weight = double.Parse(returnPositionText(s, 9));
					bexEpitope.parameters.Isoelectric_point_pH = double.Parse(returnPositionText(s, 11));
					bexEpitope.parameters.Net_Charge_At_pH_7 = double.Parse(returnPositionText(s, 12));
					if (returnPositionText(s, 13) == "Good")
					{
						bexEpitope.parameters.Innovagen_Solubility = true;
						bexEpitope.parameters.GRAVY_Score = double.Parse(s.Substring(s.LastIndexOf("\t") + 1));
					}
					else
					{
						bexEpitope.parameters.GRAVY_Score = double.Parse(returnPositionText(s, 14));
						bexEpitope.parameters.Innovagen_Solubility = true;
						bexEpitope.parameters.Modified_Sequence_Parameters = new ModifiedParameters();
						bexEpitope.parameters.Modified_Sequence_Parameters.Modified_Sequence = returnPositionText(s, 15);
						bexEpitope.parameters.Modified_Sequence_Parameters.GRAVY_Score = double.Parse(returnPositionText(s, 16));
					}
					//find the allele data
					foreach (Protein p in deserializedDatas[2].all)
					{
						if (p.proteinName == bexProtein.proteinName)
						{
							foreach (Epitope e in p.epitopes)
							{
								if (e.epitopeSeq == bexEpitope.epitopeSeq)
								{
									bexEpitope.alleles = e.alleles;
								}
							}
						}
					}

					bexEpList.Add(bexEpitope);
					bexEpitope = new Epitope();
				}
			}
			#endregion
			#region outputtable
			//bexData = 2.2, bexData2 = iedb.
			//table format: Name > Allegen > Epitope > Cores (2.2) > start > Tepitope score.
			if (!Directory.Exists(resultDir + "/table"))
				Directory.CreateDirectory(resultDir + "/table");
			//foreach (Protein s in bexData.all)
			for (int bd = 0; bd < bexData.all.Count; bd++)
			{
				if (bexData.all[bd].proteinName != null)
				{
				string path = resultDir + "/table/" + bexData.all[bd].proteinName + ".xlsx";
				//excel shet.
				/*
                Excel.Application excelApp = new Excel.Application();
                Excel.Workbook excelWorkbook1;
                Excel.Range excelCell;
                //save a new xslx file
                excelWorkbook1 = excelApp.Workbooks.Add();
                excelWorkbook1.SaveAs(path);
                excelWorkbook1.Close();
                excelWorkbook1 = excelApp.Workbooks.Open(path, 0, false, 5, "", "HelloKitty", false, Excel.XlPlatform.xlWindows, "", true, false, 0, true, false, false);

                //opens the first worksheet.
                Excel.Sheets excelSheets1 = excelWorkbook1.Worksheets;
                Excel.Worksheet excelWorksheet1 = (Excel.Worksheet)excelSheets1[1];
                excelWorksheet1.Cells.Style.HorizontalAlignment = Microsoft.Office.Interop.Excel.XlHAlign.xlHAlignCenter;
                excelWorksheet1.Cells.Style.VerticalAlignment = Microsoft.Office.Interop.Excel.XlVAlign.xlVAlignCenter;
                */

				ExcelPackage pck = new ExcelPackage(new FileInfo(path));
				ExcelWorksheet ws = pck.Workbook.Worksheets.Add(bexData.all[bd].proteinName);
				//insert data here. 
				//do headings:
				#region table 1
				ws.Cells["A1"].Value = "Proposed Name";
				ws.Cells["B1"].Value = "Allergen";
				ws.Cells["C1"].Value = "Epitope Sequence";
				ws.Cells["D1"].Value = "Binding Cores";
				ws.Cells["E1"].Value = "Position";
				ws.Cells["F1"].Value = "*0101 (7.27%)";
				ws.Cells["G1"].Value = "*0301 (12.77%)";
				ws.Cells["H1"].Value = "*0401 (6.54%)";
				ws.Cells["I1"].Value = "*0404 (3.05%)";
				ws.Cells["J1"].Value = "*0405 (11.70%)";
				ws.Cells["K1"].Value = "*0701 (11.70%)";
				ws.Cells["L1"].Value = "*0802 (0.26%)";
				ws.Cells["M1"].Value = "*1101 (5.28%)";
				ws.Cells["N1"].Value = "*1302 (4.10%)";
				ws.Cells["O1"].Value = "*1501 (11.05%)";

				/*
                excelCell = (Excel.Range)excelWorksheet1.get_Range("A1"); 
                excelCell.Value2 = "Proposed Name";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("B1");
                excelCell.Value2 = "Allergen";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("C1");
                excelCell.Value2 = "Epitope Sequence";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("D1");
                excelCell.Value2 = "Binding Cores";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("E1");
                excelCell.Value2 = "Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("F1");
                excelCell.Value2 = "*0101 (7.27%)";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("G1");
                excelCell.Value2 = "*0301 (12.77%)"; //"Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("H1");
                excelCell.Value2 = "*0401 (6.54%)";// "Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("I1");
                excelCell.Value2 = "*0404 (3.05%)"; //"Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("J1");
                excelCell.Value2 = "*0405 (11.70%)"; //"Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("K1");
                excelCell.Value2 = "*0701 (11.70%)"; //"Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("L1");
                excelCell.Value2 = "*0802 (0.26%)"; //"Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("M1");
                excelCell.Value2 = "*1101 (5.28%)"; //"Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("N1");
                excelCell.Value2 = "*1302 (4.10%)"; //"Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("O1");
                excelCell.Value2 = "*1501 (11.05%)"; //"Start position";*/

				int row, column, counter;
				row = 2; column = 1; counter = 1;
				foreach (Epitope e in bexData2.all[bd].epitopes)
				{
					//  excelCell = (Excel.Range)excelWorksheet1.get_Range("A" + row.ToString());
					//  excelCell.Value2 = bexData2.all[bd].proteinName + "_" + counter.ToString();
					ws.Cells["A" + row.ToString()].Value = bexData2.all[bd].proteinName + "_" + counter.ToString();
					counter++;
					//excelCell = (Excel.Range)excelWorksheet1.get_Range("B" + row.ToString());
					ws.Cells["B" + row.ToString()].Value = bexData2.all[bd].proteinName;
					//  excelCell = (Excel.Range)excelWorksheet1.get_Range("C" + row.ToString());
					//if (e.parameters.Modified_Sequence_Parameters != null)
					//   ws.Cells["C" + row.ToString()].Value = e.parameters.Modified_Sequence_Parameters.Modified_Sequence;
					//   else
					//      ws.Cells["C" + row.ToString()].Value = e.epitopeSeq;
					ws.Cells["C" + row.ToString()].Value = e.epitopeSeq;
					//  excelCell = (Excel.Range)excelWorksheet1.get_Range("D" + row.ToString());
					ws.Cells["D" + row.ToString()].Value = string.Join(", ",e.Binding_Core);
					//if (e.Binding_Core.Contains(","))
					//foreach(string s in e.Binding_Core)

					// else
					//    ws.Cells["D" + row.ToString()].Value = e.Binding_Core;
					// excelCell = (Excel.Range)excelWorksheet1.get_Range("E" + row.ToString());
					ws.Cells["E" + row.ToString()].Value = e.parameters.Epitope_Start.ToString() + "~" + (e.parameters.Epitope_Start+14).ToString();
					//   excelCell = (Excel.Range)excelWorksheet1.get_Range("F" + row.ToString());
					/* foreach (Allele a in e.alleles)
                    {
                        if (a.alleleName == "HLA-DRB1*01:01")
                        {
                            ws.Cells["F" + row.ToString()].Value = a.IC50;
                            if ((int)ws.Cells["F" + row.ToString()].Value < 1)
                                ws.Cells["F" + row.ToString()].Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
                            else if ((int)ws.Cells["F" + row.ToString()].Value < 2)
                                ws.Cells["F" + row.ToString()].Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
                            else
                                ws.Cells["F" + row.ToString()].Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
                        }
                    }
                    excelCell = (Excel.Range)excelWorksheet1.get_Range("G" + row.ToString());*/
					foreach (Allele a in e.alleles)
					{
						if (a.alleleName == "HLA-DRB1*01:01")
						{
							using (var c = ws.Cells["F" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 1)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
								else if ((double)c.Value < 2)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
							}
						}
						else if (a.alleleName == "HLA-DRB1*03:01")
						{
							using (var c = ws.Cells["G" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 1)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
								else if ((double)c.Value < 2)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
							}
						}
						else if (a.alleleName == "HLA-DRB1*04:01")
						{
							using (var c = ws.Cells["H" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 1)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
								else if ((double)c.Value < 2)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
							}
						}
						else if (a.alleleName == "HLA-DRB1*04:04")
						{
							using (var c = ws.Cells["I" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 1)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
								else if ((double)c.Value < 2)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
							}
						}
						else if (a.alleleName == "HLA-DRB1*04:05")
						{
							using (var c = ws.Cells["J" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 1)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
								else if ((double)c.Value < 2)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
							}
						}
						else if (a.alleleName == "HLA-DRB1*07:01")
						{
							using (var c = ws.Cells["K" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 1)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
								else if ((double)c.Value < 2)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
							}
						}
						else if (a.alleleName == "HLA-DRB1*08:02")
						{
							using (var c = ws.Cells["L" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 1)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
								else if ((double)c.Value < 2)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
							}
						}
						else if (a.alleleName == "HLA-DRB1*11:01")
						{
							using (var c = ws.Cells["M" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 1)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
								else if ((double)c.Value < 2)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
							}
						}
						else if (a.alleleName == "HLA-DRB1*13:02")
						{
							using (var c = ws.Cells["N" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 1)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
								else if ((double)c.Value < 2)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
							}
						}
						else if (a.alleleName == "HLA-DRB1*15:01")
						{
							using (var c = ws.Cells["O" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 1)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
								else if ((double)c.Value < 2)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
							}
						}
					}
					row++;
				}

				row++;
				counter = 1;
				ws.Cells["A" + row.ToString()].Value = "Proposed Name";
				ws.Cells["B" + row.ToString()].Value = "Allergen";
				ws.Cells["C" + row.ToString()].Value = "Epitope Sequence";
				ws.Cells["D" + row.ToString()].Value = "Binding Cores";
				ws.Cells["E" + row.ToString()].Value = "Position";
				ws.Cells["F" + row.ToString()].Value = "*0101 (7.27%)";
				ws.Cells["G" + row.ToString()].Value = "*0301 (12.77%)";
				ws.Cells["H" + row.ToString()].Value = "*0401 (6.54%)";
				ws.Cells["I" + row.ToString()].Value = "*0404 (3.05%)";
				ws.Cells["J" + row.ToString()].Value = "*0405 (11.70%)";
				ws.Cells["K" + row.ToString()].Value = "*0701 (11.70%)";
				ws.Cells["L" + row.ToString()].Value = "*0802 (0.26%)";
				ws.Cells["M" + row.ToString()].Value = "*1101 (5.28%)";
				ws.Cells["N" + row.ToString()].Value = "*1302 (4.10%)";
				ws.Cells["O" + row.ToString()].Value = "*1501 (11.05%)";
				/*
                excelCell = (Excel.Range)excelWorksheet1.get_Range("A" + row.ToString());
                excelCell.Value2 = "Proposed Name";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("B" + row.ToString());
                excelCell.Value2 = "Allergen";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("C" + row.ToString());
                excelCell.Value2 = "Epitope Sequence";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("D" + row.ToString());
                excelCell.Value2 = "Binding Cores";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("E" + row.ToString());
                excelCell.Value2 = "Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("F" + row.ToString());
                excelCell.Value2 = "*0101 (7.27%)";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("G" + row.ToString());
                excelCell.Value2 = "*0301 (12.77%)"; //"Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("H" + row.ToString());
                excelCell.Value2 = "*0401 (6.54%)";// "Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("I" + row.ToString());
                excelCell.Value2 = "*0404 (3.05%)"; //"Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("J" + row.ToString());
                excelCell.Value2 = "*0405 (11.70%)"; //"Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("K" + row.ToString());
                excelCell.Value2 = "*0701 (11.70%)"; //"Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("L" + row.ToString());
                excelCell.Value2 = "*0802 (0.26%)"; //"Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("M" + row.ToString());
                excelCell.Value2 = "*1101 (5.28%)"; //"Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("N" + row.ToString());
                excelCell.Value2 = "*1302 (4.10%)"; //"Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("O" + row.ToString());
                excelCell.Value2 = "*1501 (11.05%)"; //"Start position";
                if (bexData == bexData2)
                {

                }*/
				row++;
				foreach (Epitope e2 in bexData.all[bd].epitopes)
				{
					// excelCell = (Excel.Range)excelWorksheet1.get_Range("A" + row.ToString());
					ws.Cells["A" + row.ToString()].Value = bexData.all[bd].proteinName + "_" + counter.ToString();
					counter++;
					// excelCell = (Excel.Range)excelWorksheet1.get_Range("B" + row.ToString());
					ws.Cells["B" + row.ToString()].Value = bexData.all[bd].proteinName;
					// excelCell = (Excel.Range)excelWorksheet1.get_Range("C" + row.ToString());
					// if (e2.parameters.Modified_Sequence_Parameters != null)
					//     ws.Cells["C" + row.ToString()].Value = e2.parameters.Modified_Sequence_Parameters.Modified_Sequence;
					//   else
					//      ws.Cells["C" + row.ToString()].Value = e2.epitopeSeq;
					ws.Cells["C" + row.ToString()].Value = e2.epitopeSeq;
					ws.Cells["D" + row.ToString()].Value = string.Join(", ", e2.Binding_Core);
					// excelCell = (Excel.Range)excelWorksheet1.get_Range("D" + row.ToString());
					//  if (e2.Binding_Core.Contains(","))
					//     ws.Cells["D" + row.ToString()].Value = e2.Binding_Core.Replace(",", ",");
					//  else
					//      ws.Cells["D" + row.ToString()].Value = e2.Binding_Core;
					//  excelCell = (Excel.Range)excelWorksheet1.get_Range("E" + row.ToString());
					ws.Cells["E" + row.ToString()].Value = e2.parameters.Epitope_Start.ToString() + "~" + (e2.parameters.Epitope_Start + 14).ToString();

					foreach (Allele a in e2.alleles)
					{
						if (a.alleleName == "HLA-DRB10101")
						{
							using (var c = ws.Cells["F" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 50)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
								else if ((double)c.Value < 500)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
							}
						}
						else if (a.alleleName == "HLA-DRB10301")
						{
							using (var c = ws.Cells["G" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 50)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
								else if ((double)c.Value < 500)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
							}
						}
						else if (a.alleleName == "HLA-DRB10401")
						{
							using (var c = ws.Cells["H" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 50)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
								else if ((double)c.Value < 500)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
							}
						}
						else if (a.alleleName == "HLA-DRB10404")
						{
							using (var c = ws.Cells["I" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 50)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
								else if ((double)c.Value < 500)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
							}
						}
						else if (a.alleleName == "HLA-DRB10405")
						{
							using (var c = ws.Cells["J" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 50)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
								else if ((double)c.Value < 500)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
							}
						}
						else if (a.alleleName == "HLA-DRB10701")
						{
							using (var c = ws.Cells["K" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 50)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
								else if ((double)c.Value < 500)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
							}
						}
						else if (a.alleleName == "HLA-DRB10802")
						{
							using (var c = ws.Cells["L" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 50)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
								else if ((double)c.Value < 500)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
							}
						}
						else if (a.alleleName == "HLA-DRB11101")
						{
							using (var c = ws.Cells["M" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 50)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
								else if ((double)c.Value < 500)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
							}
						}
						else if (a.alleleName == "HLA-DRB11302")
						{
							using (var c = ws.Cells["N" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 50)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
								else if ((double)c.Value < 500)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
							}
						}
						else if (a.alleleName == "HLA-DRB11501")
						{
							using (var c = ws.Cells["O" + row.ToString()])
							{
								c.Value = a.IC50;
								c.Style.Fill.PatternType = ExcelFillStyle.Solid;
								if ((double)c.Value < 50)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Green);
								else if ((double)c.Value < 500)
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Yellow);
								else
									c.Style.Fill.BackgroundColor.SetColor(System.Drawing.Color.Red);
							}
						}
					}
					row++;
				}

				row++;
				counter = 1;
				//Epitope Name    Overall Score   Epitope Sequence    Binding Cores   Start Position  GRAVY Score Isoelectric Point   Soluble Sequence (GRAVY Score)
				/*excelCell = (Excel.Range)excelWorksheet1.get_Range("A" + row.ToString());
                excelCell.Value2 = "Proposed name";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("B" + row.ToString());
                excelCell.Value2 = "Soluble sequence";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("C" + row.ToString());
                excelCell.Value2 = "Epitope start";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("D" + row.ToString());
                excelCell.Value2 = "Isoelectric point";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("E" + row.ToString());
                excelCell.Value2 = "GRAVY score";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("F" + row.ToString());
                excelCell.Value2 = "Water solubility";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("G" + row.ToString());
                excelCell.Value2 = "Tepitope score";// "Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("H" + row.ToString());
                excelCell.Value2 = "NetMHCII score"; //"Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("I" + row.ToString());
                excelCell.Value2 = "Overall score";// "Start position";
                excelCell = (Excel.Range)excelWorksheet1.get_Range("J" + row.ToString());*/

				ws.Cells["A" + row.ToString()].Value = "Proposed name";
				ws.Cells["B" + row.ToString()].Value = "Soluble sequence";
				ws.Cells["C" + row.ToString()].Value = "Epitope Position";
				ws.Cells["D" + row.ToString()].Value = "Isoelectric point";
				ws.Cells["E" + row.ToString()].Value = "GRAVY score";
				ws.Cells["F" + row.ToString()].Value = "Water solubility";
				ws.Cells["G" + row.ToString()].Value = "Tepitope score";
				ws.Cells["H" + row.ToString()].Value = "NetMHCII2.2 score";
				ws.Cells["I" + row.ToString()].Value = "Modified NetMHCII2.2 score";
				ws.Cells["J" + row.ToString()].Value = "Overall score";
				ws.Cells["K" + row.ToString()].Value = "Rank";
				row++;

				foreach (Epitope e in bexData.all[bd].epitopes)
				{
					foreach (Protein p in deserializedDatas[0].all)
					{
						if (p.proteinName == bexData.all[bd].proteinName)
						{
							foreach (Epitope ep in p.epitopes)
							{
								if (ep.epitopeSeq == e.epitopeSeq)
								{
									// excelCell = (Excel.Range)excelWorksheet1.get_Range("H" + row.ToString());
									ws.Cells["H" + row.ToString()].Value = ep.average.ToString("F");
								}
							}
						}
					}
					foreach (Protein p in deserializedDatas[2].all)
					{
						if (p.proteinName == bexData.all[bd].proteinName)
						{
							foreach (Epitope ep in p.epitopes)
							{
								if (ep.epitopeSeq == e.epitopeSeq)
								{
									// excelCell = (Excel.Range)excelWorksheet1.get_Range("G" + row.ToString());
									ws.Cells["G" + row.ToString()].Value = ep.average.ToString("F");// = ep.average.ToString("F");
								}
							}
						}
					}
					//Epitope Name    Overall Score   Epitope Sequence    Binding Cores   Start Position  GRAVY Score Isoelectric Point   Soluble Sequence (GRAVY Score)
					// excelCell = (Excel.Range)excelWorksheet1.get_Range("A" + row.ToString());
					ws.Cells["A" + row.ToString()].Value = bexData.all[bd].proteinName + "_" + counter.ToString();
					//  excelCell = (Excel.Range)excelWorksheet1.get_Range("B" + row.ToString());

					//excelWorksheet1.Cells[2,row].IsRichText = true;
					if (e.parameters.Modified_Sequence_Parameters != null)
						ws.Cells["B" + row.ToString()].Value = e.parameters.Modified_Sequence_Parameters.Modified_Sequence;
					else
						ws.Cells["B" + row.ToString()].Value = e.epitopeSeq;
					//  excelCell = (Excel.Range)excelWorksheet1.get_Range("C" + row.ToString());
					ws.Cells["C" + row.ToString()].Value = e.parameters.Epitope_Start.ToString() + "~" + (e.parameters.Epitope_Start + 14).ToString();
					//  excelCell = (Excel.Range)excelWorksheet1.get_Range("D" + row.ToString())
					if (e.parameters.Modified_Sequence_Parameters != null)
						ws.Cells["D" + row.ToString()].Value = e.parameters.Modified_Sequence_Parameters.Isoelectric_point_pH;
					else
						ws.Cells["D" + row.ToString()].Value = e.parameters.Isoelectric_point_pH;
					//excelCell = (Excel.Range)excelWorksheet1.get_Range("E" + row.ToString());
					if (e.parameters.Modified_Sequence_Parameters != null)
						ws.Cells["E" + row.ToString()].Value = e.parameters.Modified_Sequence_Parameters.GRAVY_Score;
					else
						ws.Cells["E" + row.ToString()].Value = e.parameters.GRAVY_Score;
					if (e.parameters.Modified_Sequence_Parameters != null)
						ws.Cells["I" + row.ToString()].Value = e.parameters.Modified_Sequence_Parameters.Modified_Average;
					else
						ws.Cells["I" + row.ToString()].Value = "N/A";
					// excelCell = (Excel.Range)excelWorksheet1.get_Range("F" + row.ToString());
					ws.Cells["F" + row.ToString()].Value = "Good";
					//   excelCell = (Excel.Range)excelWorksheet1.get_Range("I" + row.ToString());
					ws.Cells["J" + row.ToString()].Value = e.average;// "Start position";
					// excelCell = (Excel.Range)excelWorksheet1.get_Range("J" + row.ToString());
					ws.Cells["K" + row.ToString()].Value = counter; //"Start position";

					row++;
					counter++;
				}


				// excelWorksheet1.Columns.AutoFit();
				//  excelWorksheet1.Rows.AutoFit();

				//  excelWorkbook1.Save();
				//  excelWorkbook1.Close();
				// ws.Column(1).AutoFit(0);
				ws.Cells.AutoFitColumns();
				pck.Save();
				l.writeLog("[INFO] successfully wrote to excel file: " + path);
				#endregion

				//modify that file
				FileInfo fi = new FileInfo(path);
				using (ExcelPackage package = new ExcelPackage(fi))
				{
					ExcelWorksheet worksheet = package.Workbook.Worksheets[bexData.all[bd].proteinName];
					for (int x = 1; x < 9999999; x++)
					{
						if (worksheet.Cells[x, 2].Value == null)
						{

						}
						else if (worksheet.Cells[x,2].Value.ToString() == "Soluble sequence")
						{
							x++;
							int count2 = -1;
							while (worksheet.Cells[x, 2].Value != null) 
							{
								count2 ++;
								//string seq = worksheet.Cells[x, 2].Value.ToString();
								worksheet.Cells[x, 2].Value = "";
								string pre = "";
								string post = "";
								int corestart = 99;
								int coreend = 0;
								string precore = "";
								string postcore = "";
								string core = "";

								//using (Epitope e = bexData.all[bd].epitopes[count2])
								//foreach (Epitope e in bexData.all[bd].epitopes)
								//{
									string seq = "";//e.epitopeSeq;
									List<string> cores = bexData.all[bd].epitopes[count2].Binding_Core;
									if (bexData.all[bd].epitopes[count2].parameters.Modified_Sequence_Parameters != null) //modified
									{
										seq = bexData.all[bd].epitopes[count2].parameters.Modified_Sequence_Parameters.Modified_Sequence;
										if (bexData.all[bd].epitopes[count2].epitopeSeq != seq)
										{
											if (seq.IndexOf(bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S")) == -1)
											{
												if (!seq.Contains("Unable"))
												{

												}
												pre = "!NF"; post = "!NF";
											}
											else if (seq.IndexOf(bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S")) == 0)
											{
												pre = "";
												post = seq.Substring(bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S").Length);
											}
											else if (seq.IndexOf(bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S")) + bexData.all[bd].epitopes[count2].epitopeSeq.Length == seq.Length) //nothing after
											{
												pre = seq.Substring(0, seq.IndexOf(bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S")));
												post = "";
											}
											else
											{
												pre = seq.Substring(0, seq.IndexOf(bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S")));
												post = seq.Substring(pre.Length + bexData.all[bd].epitopes[count2].epitopeSeq.Length);
											}

											foreach (string c in cores)
											{
												if (bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S").IndexOf(c) != -1)
												{
													if (bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S").IndexOf(c) < corestart)
														corestart = bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S").IndexOf(c);
													if (bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S").IndexOf(c) + 9 > coreend)
														coreend = bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S").IndexOf(c) + 9;
												}
												else
												{

												}
											}

											if (corestart < coreend)
											{
												if (corestart == 0 && coreend == bexData.all[bd].epitopes[count2].epitopeSeq.Length)
												{
													precore = "";
													postcore = "";
												}
												else if (corestart == 0)
												{
													precore = "";
													postcore = bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S").Substring(coreend);
												}
												else if (coreend == bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S").Length)
												{
													precore = bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S").Substring(0, corestart);
													postcore = "";
												}
												else
												{
													precore = bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S").Substring(0, corestart);
													postcore = bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S").Substring(coreend);
												}
												core = seq.Substring(corestart + pre.Length, coreend - corestart);
											}
											else
											{

											}
											if (core != "")
											{
												worksheet.Cells[x, 2].IsRichText = true;
												ExcelRichTextCollection rtfCollection = worksheet.Cells[x, 2].RichText;
												ExcelRichText ert;
												if (pre != "")
												{
													ert = rtfCollection.Add(pre);
													ert.Bold = true;
													ert.Italic = false;
													ert.Color = System.Drawing.Color.Red;
												}
												if (precore != "")
												{
													ert = rtfCollection.Add(precore);
													ert.Bold = false;
													ert.Italic = false;
													ert.Color = System.Drawing.Color.Black;
												}
												ert = rtfCollection.Add(core);
												ert.Bold = false;
												ert.Italic = true;
												ert.Color = System.Drawing.Color.Blue;
												if (postcore != "")
												{
													ert = rtfCollection.Add(postcore);
													ert.Bold = false;
													ert.Italic = false;
													ert.Color = System.Drawing.Color.Black;
												}
												if (post != "")
												{
													ert = rtfCollection.Add(post);
													ert.Bold = true;
													ert.Italic = false;
													ert.Color = System.Drawing.Color.Red;
												}
											}
											else
											{

											}
										}
									}
									else if (bexData.all[bd].epitopes[count2].parameters.Modified_Sequence_Parameters == null)
									{
										//if (bexData.all[bd].epitopes[count2].epitopeSeq.Replace("C", "S") == seq)
										//{
											foreach (string c in cores)
											{
												if (bexData.all[bd].epitopes[count2].epitopeSeq.IndexOf(c) != -1)
												{
													if (bexData.all[bd].epitopes[count2].epitopeSeq.IndexOf(c) < corestart)
														corestart = bexData.all[bd].epitopes[count2].epitopeSeq.IndexOf(c);
													if (bexData.all[bd].epitopes[count2].epitopeSeq.IndexOf(c) + 9 > coreend)
														coreend = bexData.all[bd].epitopes[count2].epitopeSeq.IndexOf(c) + 9;
												}
												else
												{

												}
											}

											if (corestart < coreend)
											{
												if (corestart == 0 && coreend == bexData.all[bd].epitopes[count2].epitopeSeq.Length)
												{
													precore = "";
													postcore = "";
												}
												else if (corestart == 0)
												{
													precore = "";
													postcore = bexData.all[bd].epitopes[count2].epitopeSeq.Substring(coreend);
												}
												else if (coreend == bexData.all[bd].epitopes[count2].epitopeSeq.Length)
												{
													precore = bexData.all[bd].epitopes[count2].epitopeSeq.Substring(0, corestart);
													postcore = "";
												}
												else
												{
													precore = bexData.all[bd].epitopes[count2].epitopeSeq.Substring(0, corestart);
													postcore = bexData.all[bd].epitopes[count2].epitopeSeq.Substring(coreend);
												}
											core = bexData.all[bd].epitopes[count2].epitopeSeq.Substring(corestart + pre.Length, coreend - corestart - pre.Length );
											}

											if (core != "")
											{
												worksheet.Cells[x, 2].IsRichText = true;
												ExcelRichTextCollection rtfCollection = worksheet.Cells[x, 2].RichText;
												ExcelRichText ert;
												if (precore != "")
												{
													ert = rtfCollection.Add(precore);
													ert.Italic = false;
													ert.Color = System.Drawing.Color.Black;
												}
												ert = rtfCollection.Add(core);
												ert.Italic = true;
												ert.Color = System.Drawing.Color.Blue;
												if (postcore != "")
												{
													ert = rtfCollection.Add(postcore);
													ert.Italic = false;
													ert.Color = System.Drawing.Color.Black;
												}
											}
									//  }
									}
							    //end foreach loop.
								//break;
								x++;
							} 
							break;
						}
					}
					package.Save();
				}
				}
			}
			#endregion


		}

		public void compileProteinseq()
		{
			foreach (string s in Directory.GetFiles(resultDir)) //read in the 3 json files.
			{
				if (s.Contains(".json"))
				{
					if (s.Contains("netMHCII-2.2_results.json"))
					{
						using (StreamReader f = File.OpenText(s))
						{
							JsonSerializer serializer = new JsonSerializer();
							deserializedDatas[0] = (rawData)serializer.Deserialize(f, typeof(rawData));
						}
					}
				}
			}

			List<string> output = new List<string>();
			output.Add("Allergen Name\tAllergen Accession\tAllergen Sequence");
			foreach (Protein p in deserializedDatas[0].all)
			{
				output.Add(p.proteinName + "\t" + p.proteinAccession + "\t" + p.proteinSeq);                
			}
			l.write(output, path: resultDir + "/seq.txt");
		}
	}
}
