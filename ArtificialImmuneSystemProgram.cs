using System;
using System.Collections.Generic;
using System.Collections; // para BitArray

namespace ArtificialImmuneSystem
{
  class ArtificialImmuneSystemProgram
  {
    static Random random;

    static void Main(string[] args)
    {
      Console.WriteLine("\nInicia la simulacion del Sistema Inmune Artificial\n");

      random = new Random(1);

      int numPatternBits = 12;
      int numAntibodyBits = 4;
      int numLymphocytes = 3;
      int stimulationThreshold = 3;

      Console.WriteLine("Cargando conjunto de auto-antígeno (patrones históricos 'normales')");
      List<BitArray> selfSet = LoadSelfSet(null);
      ShowSelfSet(selfSet);

      Console.WriteLine("\nCrear conjunto de linfocitos usando seleccion negativa" +
        " y deteccion r-chunks");
      List<Lymphocyte> lymphocyteSet = CreateLymphocyteSet(selfSet, numAntibodyBits,
        numLymphocytes);
      ShowLymphocyteSet(lymphocyteSet);

      Console.WriteLine("\nComience la simulacion de deteccion de intrusos AIS\n");

      int time = 0;
      int maxTime = 6;

      while (time < maxTime)
      {
        Console.WriteLine("===============================================");
        //Console.WriteLine("time = "+ time);
        BitArray incoming = RandomBitArray(numPatternBits);
        Console.WriteLine("Patron entrante = " + BitArrayAsString(incoming) + "\n");

        for (int i = 0; i < lymphocyteSet.Count; ++i)
        {
          if (lymphocyteSet[i].Detects(incoming) == true)
          {
            Console.WriteLine("Patron entrante detectado por linfocito " + i);
            ++lymphocyteSet[i].stimulation;
            if (lymphocyteSet[i].stimulation >= stimulationThreshold)
              Console.WriteLine("linfocito " + i + " estimulado!" +
                " Comprobar intrusion entrante posible!");
            else
              Console.WriteLine("linfocito " + i + " no estimulado sobre el umbral");
          }
          else
            Console.WriteLine("Patron de entrada no detectado por el linfocito " + i);
        }
        ++time;
        Console.WriteLine("===============================================");
      }  // while

      Console.WriteLine("\nFIN AIS IDS\n");
      Console.ReadLine();
    }  // Main

    public static List<BitArray> LoadSelfSet(string dataSource)
    {
      List<BitArray> result = new List<BitArray>();
      bool[] self0 = new bool[] { true, false, false, true, false, true, true, false, true, false, false, true };
      bool[] self1 = new bool[] { true, true, false, false, true, false, true, false, true, true, false, false };
      bool[] self2 = new bool[] { true, false, true, true, false, false, true, true, false, true, false, true };
      bool[] self3 = new bool[] { false, false, true, true, false, true, false, true, true, false, true, true };
      bool[] self4 = new bool[] { false, true, false, true, false, true, false, false, true, true, false, true };
      bool[] self5 = new bool[] { false, false, true, false, true, false, true, false, false, true, false, false };

      result.Add(new BitArray(self0));
      result.Add(new BitArray(self1));
      result.Add(new BitArray(self2));
      result.Add(new BitArray(self3));
      result.Add(new BitArray(self4));
      result.Add(new BitArray(self5));

      return result;
    }

    public static void ShowSelfSet(List<BitArray> selfSet)
    {
      for (int i = 0; i < selfSet.Count; ++i)
      {
        Console.WriteLine(i + ": " + BitArrayAsString(selfSet[i]));
      }
    }

    public static string BitArrayAsString(BitArray ba)
    {
      string s = "";
      for (int i = 0; i < ba.Length; ++i)
        s += (ba[i] == true) ? "1 " : "0 ";
      return s;
    }

    public static List<Lymphocyte> CreateLymphocyteSet(List<BitArray> selfSet, int numAntibodyBits, int numLymphocytes)
    {
	// crea una lista de linfocitos que no detecta ningún patrón
      List<Lymphocyte> result = new List<Lymphocyte>();
      Dictionary<int, bool> contents = new Dictionary<int, bool>();  

      while (result.Count < numLymphocytes)
      {
        BitArray antibody = RandomBitArray(numAntibodyBits);  // random anticuerpo
        Lymphocyte lymphocyte = new Lymphocyte(antibody);     // random linfocito
        int hash = lymphocyte.GetHashCode();                  // asume la longitud del anticuerpo <= 32 bits

        if (DetectsAny(selfSet, lymphocyte) == false && contents.ContainsKey(hash) == false)
        {
          result.Add(lymphocyte);
          contents.Add(hash, true);
        }
      }
      return result;
    }

    private static bool DetectsAny(List<BitArray> selfSet, Lymphocyte lymphocyte)  // helper
    {

      for (int i = 0; i < selfSet.Count; ++i)
      {
        if (lymphocyte.Detects(selfSet[i]) == true)
          return true;
      }
      return false;
    }

    public static void ShowLymphocyteSet(List<Lymphocyte> lymphocyteySet)
    {
      for (int i = 0; i < lymphocyteySet.Count; ++i)
      {
        Console.WriteLine(i + ": " + lymphocyteySet[i].ToString());
      }
    }

    public static BitArray RandomBitArray(int numBits)
    {
      bool[] bools = new bool[numBits];
      for (int i = 0; i < numBits; ++i)
      {
        int b = random.Next(0, 2);
        bools[i] = (b == 0) ? false : true;
      }
      return new BitArray(bools);
    }

  }  

  public class Lymphocyte
  {
    public BitArray antibody;  // Detector
    public int[] searchTable;  // Para detección rápida
    public int age;            // Determina su destrucción
    public int stimulation;    // Activador

    public Lymphocyte(BitArray antibody)
    {
      this.antibody = new BitArray(antibody); // asume len >= 2
      this.searchTable = BuildTable();        
      this.age = 0;
      this.stimulation = 0;
    }

    private int[] BuildTable()
    {
      int[] result = new int[antibody.Length];
      int pos = 2;
      int cnd = 0;
      result[0] = -1;
      result[1] = 0;
      while (pos < antibody.Length)
      {
        if (antibody[pos - 1] == antibody[cnd])
        {
          ++cnd; result[pos] = cnd; ++pos;
        }
        else if (cnd > 0)
          cnd = result[cnd];
        else
        {
          result[pos] = 0; ++pos;
        }
      }
      return result;
    }

    public bool Detects(BitArray pattern)
    {
      // does the this.antibody detector detect pattern?
      // Knuth-Morris-Pratt algorithm aka r-chunks
      
      int m = 0;
      int i = 0;
      while (m + i < pattern.Length)
      {
        if (this.antibody[i] == pattern[m + i])
        {
          if (i == antibody.Length - 1)
            return true;
          ++i;
        }
        else
        {
          m = m + i - this.searchTable[i];
          if (searchTable[i] > -1)
            i = searchTable[i];
          else
            i = 0;
        }
      }
      return false;  // not found
    }

    public override int GetHashCode()
    {
      int[] singleInt = new int[1];
      antibody.CopyTo(singleInt, 0);
      return singleInt[0];
    }

    public override string ToString()
    {
      string s = "antibody = ";
      for (int i = 0; i < antibody.Length; ++i)
        s += (antibody[i] == true) ? "1 " : "0 ";
      s += " age = " + age;
      s += "  stimulation = " + stimulation;
      return s;
    }
   
  }  

}  
