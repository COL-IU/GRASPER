/**
>HEADER
    Copyright (c) 2014/2015 Heewook Lee heewlee@indiana.edu

    This file is part of the GRASPER suite.

    GRASPER is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GRASPER is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GRASPER.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/

import java.io.*;

public class Constants{

    public static int MEDIAN = 422;
    
    public static int MAD = 59;

    public static int STD = Constants.MAD * 3 / 2;

    public static int NSTD = 5; //for clustering

    public static int DISCORDANT_NSTD = 5; //for filtering

    public static int MAX_D_CORDANT = Constants.MEDIAN + (Constants.DISCORDANT_NSTD * Constants.STD);

    /* this is the MAX ENCOMPASSING SIZE OF PAIRED_END READS */
    public static int MAX_D = Constants.MEDIAN + (Constants.NSTD * Constants.STD);
   
    /* this is the read length of reads*/
    public static int READLEN = 90;

    public static int MIN_MAPPINGLEN = 60;

    public static int MAX_D_MIDPOINTS = Constants.MAX_D - Constants.READLEN;

    public static int MAX_D_MIDPOINTS_CORDANT = Constants.MAX_D_CORDANT - Constants.READLEN;

    /* this is the dm in clustering algorithm but it's readlen shorted as this is adjusted for inter-midpoint distance */
    public static int DM = Constants.MAX_D - (2*Constants.READLEN) + Constants.READLEN/2; //ADDED the last term (Constants.READLEN/2) to allow overhang 
    
    /* SHOULD NOT BE USED as BWA PUTS 0 for repetitive mapping result --  MAP QUALITY for mapping to use */
    public static int MIN_MAP_QUAL = 0;
    
    /* MIN DISTANCE */
    public static int MIN_D_CORDANT = Constants.MAX_D_CORDANT - (Constants.DISCORDANT_NSTD * 2 * Constants.STD); // *2since we are subtracting from MAX_D

    public static int MIN_D = Constants.MAX_D - (Constants.NSTD * 2 * Constants.STD); // *2since we are subtracting from MAX_D

    public static int MIN_D_MIDPOINTS = Constants.MIN_D - Constants.READLEN;

    public static int MIN_D_MIDPOINTS_CORDANT = Constants.MIN_D_CORDANT - Constants.READLEN;

    public static int MIN_COVERAGE = 5;
    
    public static int MIN_NUM_READS_PER_CLUSTER = 5;

    /* MAX DISTANCE == DM */
    //public static int MAX_D = Constants.DM;// Constants.DM + (Constants.NSTD * Constants.STD);
        
    public static boolean DEBUG = false;//true;
    
    public static boolean DEBUG2 = false;
    
    public static boolean DEBUG3 = false;
    
    public static boolean DEBUG4 = false;
    
    public static boolean DYNAMIC_READLEN = true;

    /* GENOMELEN */
    public static int GENOMELEN = 3284156;//4639675;
    
    public static int DEPTH_WIN_SIZE = 100;
    
    public static int MIN_DEPTH = 1;

    public static double NO_COVERAGE_RATIO = 0.9;

    public static int CLUSTER_OVERLAP = -50;
    
    public static int MAX_CLUSTER_GAP = 50;

    public static int MIN_INVERSION_SIZE = 300;
    
    public static int MAX_INVERSION_SIZE = 1000000;

    public static int MIN_TRANSPOSITION_TARGET = 200;
    public static int MAX_TRANSPOSITION_TARGET = 1000000;

    public static int MIN_DELETION_SIZE = 100;
    public static int MIN_DELETION_SPANNING_SIZE_INCLUDING_CLUSTERS = Constants.MIN_DELETION_SIZE + MIN_D*2;
    
    public static int MIN_TANDUP_SIZE = 1000;
    public static int MAX_TANDUP_SIZE = 100000;

    public static String PROJNAME = "PROJNAME";

    public static String MEDMADFILE = "";

    public static String REFSEQ = "";

    public static void printConstants(){
	System.err.println("PROJNAME:\t" + PROJNAME);
	System.err.println("REFSEQ:\t" + REFSEQ);
	System.err.println("MEDMADFILE:\t" + MEDMADFILE);
	System.err.println("MEDIAN:\t" + MEDIAN);
	System.err.println("MAD:\t" + MAD);
	System.err.println("STD:\t" + STD);
	System.err.println("NSTD:\t" + NSTD);
	System.err.println("DISCORDANT_NSTD:\t" + DISCORDANT_NSTD);
	System.err.println("MAX_D_CORDANT:\t" + MAX_D_CORDANT);
	System.err.println("MAX_D:\t" + MAX_D);
	System.err.println("READLEN:\t" + READLEN);
	System.err.println("MAX_D_MIDPOINTS:\t" + MAX_D_MIDPOINTS);
	System.err.println("MAX_D_MIDPOINTS_CORDANT:\t" + MAX_D_MIDPOINTS_CORDANT);
	System.err.println("DM:\t" + DM);
	System.err.println("MIN_MAP_QUAL:\t" + MIN_MAP_QUAL);
	System.err.println("MIN_D_CORDANT:\t" + MIN_D_CORDANT);
	System.err.println("MIN_D:\t" + MIN_D);
	System.err.println("MIN_D_MIDPOINTS_CORDANT:\t" + MIN_D_MIDPOINTS_CORDANT);
	System.err.println("MIN_COVERAGE:\t" + MIN_COVERAGE);
	System.err.println("MIN_NUM_READS_PER_CLUSTER:\t" + MIN_NUM_READS_PER_CLUSTER );
	System.err.println("DYNAMIC_READLEN:\t" + DYNAMIC_READLEN );
	System.err.println("GENOMELEN:\t" + GENOMELEN);
	System.err.println("DEPTH_WIN_SIZE:\t" + DEPTH_WIN_SIZE);
	System.err.println("MIN_DEPTH:\t" + MIN_DEPTH);
	System.err.println("NO_COVERAGE_RATIO:\t" + NO_COVERAGE_RATIO);
	System.err.println("CLUSTER_OVERLAP:\t" + CLUSTER_OVERLAP);
	System.err.println("MAX_CLUSTER_GAP:\t" + MAX_CLUSTER_GAP);
	System.err.println("MIN_TRANSPOSITION_TARGET:\t" + MIN_TRANSPOSITION_TARGET);
	System.err.println("MAX_TRANSPOSITION_TARGET:\t" + MAX_TRANSPOSITION_TARGET);
	System.err.println("MIN_DELETION_SIZE:\t" + MIN_DELETION_SIZE);
	System.err.println("MIN_DELETION_SPANNING_SIZE_INCLUDING_CLUSTERS:\t" + MIN_DELETION_SPANNING_SIZE_INCLUDING_CLUSTERS);
	System.err.println("MIN_TANDUP_SIZE:\t" + MIN_TANDUP_SIZE);
	System.err.println("MAX_TANDUP_SIZE:\t" + MAX_TANDUP_SIZE);
    }

    
    public static int compute_MAX_D_MIDPOINTS(Read r){
	if(Constants.DYNAMIC_READLEN)
	    return Constants.MAX_D - r.getLength();
	else
	    return Constants.MAX_D_MIDPOINTS;
    }

    public static int compute_MIN_D_MIDPOINTS(Read r){
	if(Constants.DYNAMIC_READLEN){
	    int tmp = Constants.MIN_D - r.getLength();
	    if(tmp < 0)
		return 0;
	    return tmp;
	}else
	    return Constants.MIN_D_MIDPOINTS;
    }

    public static int compute_MAX_D_MIDPOINTS_CORDANT(Read r){
	if(Constants.DYNAMIC_READLEN)
	    return Constants.MAX_D_CORDANT - r.getLength();
	else
	    return Constants.MAX_D_MIDPOINTS_CORDANT;
    }

    public static int compute_MIN_D_MIDPOINTS_CORDANT(Read r){
	if(Constants.DYNAMIC_READLEN){
	    int tmp = Constants.MIN_D_CORDANT - r.getLength();
	    if(tmp < 0)
		return 0;
	    return tmp;
	}else
	    return Constants.MIN_D_MIDPOINTS_CORDANT;
    }

    public static void loadConstants(int med, int mad, int nstd, int genomeLen){
    
	Constants.MEDIAN = med;
	Constants.MAD = mad;
	Constants.STD = Constants.MAD * 3 / 2; // assuming normal.
	Constants.NSTD = nstd;
	//Constants.GENOMELEN = genomeLen;
	Constants.MAX_D = Constants.MEDIAN + (Constants.NSTD * Constants.STD);
	Constants.READLEN = 90;
	Constants.MAX_D_MIDPOINTS = Constants.MAX_D - Constants.READLEN;
	Constants.DM = Constants.MAX_D - ( 2 * Constants.READLEN ); //may be possible to dynamically set this parameter
	Constants.MIN_D = Constants.MEDIAN - (Constants.NSTD * Constants.STD);
	if(Constants.MIN_D < Constants.MIN_MAPPINGLEN)
	    Constants.MIN_D = Constants.MIN_MAPPINGLEN;
	Constants.MIN_D_MIDPOINTS = Constants.MIN_D - Constants.READLEN;
	if(Constants.MIN_D_MIDPOINTS < 0)
	    Constants.MIN_D_MIDPOINTS = 0;
	Constants.MIN_COVERAGE = 5;
    }

    public static void loadConstants(String medmadFile, String configFile){
	int med = 0;
	int mad = 0;
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(medmadFile));
	    String[] tokens = br.readLine().split("\\t");
	    br.close();
	    med = (int) Double.parseDouble(tokens[0]);
	    mad = (int) Double.parseDouble(tokens[1]);
	}catch(FileNotFoundException fnfe){
	    System.err.println("medMAD file : " + medmadFile + "\tNOT FOUND!");
	    fnfe.printStackTrace();
	    System.exit(30);//using 30 for medMAD file not found.
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	Constants.loadConstants(med, mad, configFile);
    }
    
    public static void loadMedMAD(String medmadFile){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(medmadFile));
	    String[] tokens = br.readLine().split("\\t");
	    br.close();
	    Constants.MEDIAN = (int) Double.parseDouble(tokens[0]);
	    Constants.MAD = (int) Double.parseDouble(tokens[1]);
	}catch(FileNotFoundException fnfe){
	    System.err.println("medMAD file : " + medmadFile + "\tNOT FOUND!");
	    fnfe.printStackTrace();
	    System.exit(30);//using 30 for medMAD file not found.
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    public static void loadConstants(int med, int mad, String configFile){
	Constants.loadConstants(configFile, false);
	Constants.MEDIAN = med;
	Constants.MAD = mad;
	Constants.updateConstants();
	Constants.printConstants();
    }
    

    public static void loadConstants(String configFile, boolean suppressmsg){
	Constants.loadConstants(configFile, suppressmsg, true);
    }
    
    public static void loadConstants(String configFile, boolean suppressmsg, boolean loadMedMAD){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(configFile));
	    String curline = "";
	    while((curline=br.readLine())!=null){
		curline=curline.trim();
		if(!curline.equals("") && !curline.startsWith("#"))
		    Constants.processConfigLine(curline);
	    }
	    br.close();
	    if(loadMedMAD){
		if(Constants.MEDMADFILE.equals("")){
		    Constants.MEDMADFILE = Constants.PROJNAME + ".medMAD_GRASPER";
		    Constants.loadMedMAD(Constants.MEDMADFILE);
		    Constants.updateConstants();
		}
	    }else{
		Constants.updateConstants();
	    }
	    
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	if(!suppressmsg)
	    Constants.printConstants();

    }
    
    private static void loadGenomeLength(){
	if(!Constants.REFSEQ.equals("")){
	    File lenFile = new File(Constants.REFSEQ + ".len");
	    if(lenFile.exists()){
		BufferedReader br = null;
		try{
		    br = new BufferedReader(new FileReader(lenFile));
		    String curline = br.readLine();
		    if(curline!=null)
			Constants.GENOMELEN = Integer.parseInt(curline);
		    else{
			br.close();
			throw new Exception(Constants.REFSEQ + ".len : can't find genome length.\nSystem exiting..." );
		    }			
		}catch(IOException ioe){
		    System.err.println("System exiting...");
		    ioe.printStackTrace();
		    System.exit(-1);
		}catch(Exception e){
		    e.printStackTrace();
		    System.exit(-1);
		}
	    }
	}else
	    System.err.println("Can't find reference sequence.\nSystem exiting...");
    }


    public static void updateConstants(){
	Constants.STD = Constants.MAD * 3 / 2; // assuming normal.
	Constants.MAX_D = Constants.MEDIAN + (Constants.NSTD * Constants.STD);
	Constants.MAX_D_MIDPOINTS = Constants.MAX_D - Constants.READLEN;
	Constants.MAX_D_CORDANT = Constants.MEDIAN + (Constants.DISCORDANT_NSTD * Constants.STD);
	Constants.MAX_D_MIDPOINTS_CORDANT = Constants.MAX_D_CORDANT - Constants.READLEN;
	Constants.DM = Constants.MAX_D_MIDPOINTS;//Constants.MAX_D - ( 2 * Constants.READLEN ); //may be possible to dynamically set this parameter
	Constants.MIN_D = Constants.MEDIAN - (Constants.NSTD * Constants.STD);//Constants.MAX_D - (2 * Constants.NSTD * 2 * Constants.STD);
	if(Constants.MIN_D < Constants.MIN_MAPPINGLEN)
	    Constants.MIN_D = Constants.MIN_MAPPINGLEN;
	Constants.MIN_D_MIDPOINTS = Constants.MIN_D - Constants.READLEN;
	if(Constants.MIN_D_MIDPOINTS < 0)
	    Constants.MIN_D_MIDPOINTS = 0;
	Constants.MIN_D_CORDANT = Constants.MEDIAN - (Constants.DISCORDANT_NSTD * Constants.STD); 
	if(Constants.MIN_D_CORDANT < Constants.MIN_MAPPINGLEN )
	    Constants.MIN_D_CORDANT = Constants.MIN_MAPPINGLEN;
	Constants.MIN_D_MIDPOINTS_CORDANT = Constants.MIN_D_CORDANT - Constants.READLEN;
	if(Constants.MIN_D_MIDPOINTS_CORDANT < 0)
	    Constants.MIN_D_MIDPOINTS_CORDANT = 0;
	Constants.MIN_DELETION_SPANNING_SIZE_INCLUDING_CLUSTERS = Constants.MIN_DELETION_SIZE + Constants.MIN_D*2;
    }
    
    public static void processConfigLine(String line){
	if(line.startsWith("MEDIAN="))
	    Constants.MEDIAN = (int) Double.parseDouble(line.substring(line.indexOf("=")+1).trim());
	else if(line.startsWith("MAD="))
	    Constants.MAD = (int) Double.parseDouble(line.substring(line.indexOf("=")+1).trim());
	//else if(line.startsWith("STD="))
	//   Constants.STD = Integer.parseInt(line.substring(line.indexOf("=")+1).trim());
	else if(line.startsWith("NSTD="))
	    Constants.NSTD = Integer.parseInt(line.substring(line.indexOf("=")+1).trim());
	else if(line.startsWith("DISCORDANT_NSTD="))
	    Constants.DISCORDANT_NSTD = Integer.parseInt(line.substring(line.indexOf("=")+1).trim());
	else if(line.startsWith("READLEN="))
	    Constants.READLEN = Integer.parseInt(line.substring(line.indexOf("=")+1).trim());
	else if(line.startsWith("MIN_COVERAGE="))
	    Constants.MIN_COVERAGE = Integer.parseInt(line.substring(line.indexOf("=")+1).trim());
	//else if(line.startsWith("GENOMELEN="))
	//    Constants.GENOMELEN = Integer.parseInt(line.substring(line.indexOf("=")+1).trim());
	else if(line.startsWith("DYNAMIC_READLEN="))
	    Constants.DYNAMIC_READLEN = (line.substring(line.indexOf("=")+1).trim().equals("true") ? true : false);
	else if(line.startsWith("MIN_NUM_READS_PER_CLUSTER="))
	    Constants.MIN_NUM_READS_PER_CLUSTER = Integer.parseInt(line.substring(line.indexOf("=")+1).trim());
	else if(line.startsWith("MIN_TANDUP_SIZE="))
	    Constants.MIN_TANDUP_SIZE = Integer.parseInt(line.substring(line.indexOf("=")+1).trim());
	else if(line.startsWith("MAX_TANDUP_SIZE="))
	    Constants.MAX_TANDUP_SIZE = Integer.parseInt(line.substring(line.indexOf("=")+1).trim());
	else if(line.startsWith("MIN_DELETION_SIZE="))
	    Constants.MIN_DELETION_SIZE = Integer.parseInt(line.substring(line.indexOf("=")+1).trim());
	else if(line.startsWith("MIN_TRANSPOSITION_TARGET="))
	    Constants.MIN_TRANSPOSITION_TARGET = Integer.parseInt(line.substring(line.indexOf("=")+1).trim());
	else if(line.startsWith("MAX_TRANSPOSITION_TARGET="))
	    Constants.MAX_TRANSPOSITION_TARGET = Integer.parseInt(line.substring(line.indexOf("=")+1).trim());
	else if(line.startsWith("PROJNAME"))
	    Constants.PROJNAME = line.substring(line.indexOf("=")+1).trim();
	else if(line.startsWith("MEDMADFILE")) // optional (only needed if using user defined medMAD values
	    Constants.MEDMADFILE = line.substring(line.indexOf("=")+1).trim();
	else if(line.startsWith("REFSEQ")){
	    Constants.REFSEQ = line.substring(line.indexOf("=")+1).trim();
	    Constants.loadGenomeLength();
	}else
	    System.err.println("Unused config field (Skipping):\t" + line);
    }
}
