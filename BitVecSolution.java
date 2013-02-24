

import java.lang.*;
import java.util.*;
import java.io.*;

class BVNode {

  int subTreeSize;
  int parent;
  boolean match;

  BVNode() {
    subTreeSize= 0;
    parent= -2;
    match=false;
  }//Constructor

  void treeChange (int newVal){subTreeSize += newVal;}
  void parentChange (int newParent){parent = newParent;}
  void matchChange (boolean newMatch){match = newMatch;}

  int getTreeSize() {return subTreeSize;}
  int getParent() {return parent;}
  boolean getMatch() {return match;}
}//BVNode


public class BitVecSolution {

static final int BITVECTORS = 10000;
static final int BITVECLENGTH = 10000;

public static void readVectors(long[][] bitVecLongs, final int NUMOFLONGS)
//Reads the file
  {
    int i=0, j=0;
      try{
    FileInputStream fstream = new FileInputStream("bitvectors-genes.data");
    DataInputStream in = new DataInputStream(fstream);
        BufferedReader br = new BufferedReader(new InputStreamReader(in));
    String strLine;
    while ((strLine = br.readLine()) != null)   
      {
      // Read in file and assign to BitVecLongs
	for (j=0; j<NUMOFLONGS-1; j++)
	  {
	    bitVecLongs[i][j] = Long.parseLong(strLine.substring((j*64)+1, (j+1)*64), 2);
	    if (strLine.charAt(j*64)=='1')
	      {
		bitVecLongs[i][j] ^= Long.MIN_VALUE; //adds a leading '1'
	      }
	  }
        bitVecLongs[i][NUMOFLONGS-1] = Long.parseLong(strLine.substring(((NUMOFLONGS-1)*64)+1), 2);
        if (strLine.charAt((NUMOFLONGS-1)*64)=='1')
	  {
	    bitVecLongs[i][NUMOFLONGS-1] ^= Long.MIN_VALUE; //still use the leading bit in case BITVECLENGTH is a multiple of 64
	  }
        i++;
      }
    in.close();
    }catch (Exception e){System.err.println("Error: " + e.getMessage());}
    System.out.println("File has been read.");
    System.out.println(i+" lines read");
}//readVectors

public static void writeVectors(BVNode[] bitNode)
//Writes to a file
  {
  try{
  // Create file 
  FileWriter fstream = new FileWriter("bitvectors-parents.data");
  BufferedWriter out = new BufferedWriter(fstream);
    for (int i=0; i<BITVECTORS; i++)
    {
      out.write(Integer.toString(bitNode[i].getParent()));
      out.newLine();
    }
  //Close the output stream
  out.close();
  }catch (Exception e){//Catch exception if any
  System.err.println("Error: " + e.getMessage());
  }
}//writeVectors

public static int findMatch(long[][] bitVecLongs, BVNode[] bitNode, int bitIndex, final int NUMOFLONGS)
//Finds all parent-child relationships and recursively builds a family tree around an arbitrary progenitor
{
  int j,k;
  int hamDist=0;

  bitNode[bitIndex].treeChange(1);
    for (j=0; j<BITVECTORS; j++)
    {
      if (!bitNode[j].getMatch())
      {
        hamDist=0;
	hamDist += Long.bitCount(bitVecLongs[bitIndex][0] ^ bitVecLongs[j][0]);
        if (hamDist < 31)
        {
	  hamDist += Long.bitCount(bitVecLongs[bitIndex][1] ^ bitVecLongs[j][1]);
          if (hamDist < 51)
	  {
	    hamDist += Long.bitCount(bitVecLongs[bitIndex][2] ^ bitVecLongs[j][2]);
            if (hamDist < 68)
	    {
	      hamDist += Long.bitCount(bitVecLongs[bitIndex][3] ^ bitVecLongs[j][3]);
              if (hamDist < 85)
	      {
	        hamDist += Long.bitCount(bitVecLongs[bitIndex][4] ^ bitVecLongs[j][4]);
                if (hamDist < 102)
	        {
	          for (k=5; k<12; k++){hamDist += Long.bitCount(bitVecLongs[bitIndex][k] ^ bitVecLongs[j][k]);}
	          if (hamDist < 210)
		  {
	          for (k=12; k<NUMOFLONGS; k++){hamDist += Long.bitCount(bitVecLongs[bitIndex][k] ^ bitVecLongs[j][k]);}
	          if (hamDist < 2300)//2200 for a distance of 10000x10000, 130 for 500x500
		    {
	              bitNode[j].parentChange(bitIndex);
	              bitNode[j].matchChange(true);
	              bitNode[bitIndex].treeChange(findMatch(bitVecLongs, bitNode, j, NUMOFLONGS));
		    }
		  }
		}
	      }
	    }
          }
        }
      }
    }
  return(bitNode[bitIndex].getTreeSize());
}//findMatch


public static int findRealProg(BVNode[] bitNode, int cutoff)
//Finds the BitVector best suited to be the progenitor and returns its index
{
  int i;
  int progIndex=0;
  int smallest = BITVECTORS+1;

  for (i=0; i<BITVECTORS; i++)
  {
    if ((bitNode[i].getTreeSize() >= cutoff) && (bitNode[i].getTreeSize() < smallest))
    {
      progIndex = i;
      smallest = bitNode[i].getTreeSize();
    }
    bitNode[i].treeChange(-(bitNode[i].getTreeSize()+1));
  }
  return progIndex;
}//findRealProg


public static void fixTree(BVNode[] bitNode, int nodeIndex, int newParent)
//Reverses parent-child relationships so that the proper progenitor is correctly represented
{
  if(bitNode[nodeIndex].getParent()!=-1)
    {fixTree(bitNode, bitNode[nodeIndex].getParent(), nodeIndex);}
  bitNode[nodeIndex].parentChange(newParent);
}//fixTree


public static void main(String[] args) 
{
  System.out.println("Start.");
  long startTime = System.currentTimeMillis();

  int i, j, check;
  int realProg;

if (BITVECLENGTH % 64 == 0) {i = BITVECLENGTH/64;}
  else {i = (BITVECLENGTH/64)+1;}

  final int NUMOFLONGS = i;
  boolean[][] bitVector = new boolean[BITVECTORS][BITVECLENGTH];
  long [][] bitVecLongs = new long[BITVECTORS][NUMOFLONGS];
  readVectors(bitVecLongs, NUMOFLONGS);

  //Time Check
  long elapsedTime = System.currentTimeMillis()-startTime;
  float elapsedTimeSec = elapsedTime/1000F;
  System.out.printf("%.2f Seconds to Complete.%n",elapsedTimeSec);

  BVNode[] bitNode = new BVNode[BITVECTORS];
  for (i=0; i<BITVECTORS; i++){bitNode[i] = new BVNode();}

  for (i=0; i<BITVECTORS; i++){
    //This condition should only be true once, or else a relationship has been missed
    if (!bitNode[i].getMatch()){
      bitNode[i].parentChange(-1);
      bitNode[i].matchChange(true);
      System.out.println("Building Tree.");
      j=(findMatch(bitVecLongs, bitNode, i, NUMOFLONGS))/2;
      realProg = findRealProg(bitNode, j);
      fixTree(bitNode, realProg, -1);
    }
  }

  //Write to file
  writeVectors(bitNode);

  System.out.println("Done.");

  //Time Check
  elapsedTime = System.currentTimeMillis()-startTime;
  elapsedTimeSec = elapsedTime/1000F;
  System.out.printf("%.2f Seconds to Complete.%n",elapsedTimeSec);

}//main

}//hello

