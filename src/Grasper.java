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

import java.util.*;

public class Grasper{

    public static void main(String[] args) throws Exception{
	
	if(args.length > 0){
	    String command = args[0];
	    
	    if(command.equals("depth")){
		if(args.length == 4){
		    Constants.loadConstants(args[3], false, false);
		    ClusterApp obj = new ClusterApp(new Thread(Constants.GENOMELEN, args[1]));
		    obj.run(args[2], true);
		}else{
		    System.err.print("USAGE:\t");
		    Grasper.printDepthUsage();
		}
	    }else if(command.equals("sv")){
		if(args.length == 5){
		    Constants.loadConstants(args[3], false);
		    ClusterApp obj = new ClusterApp(new Thread(Constants.GENOMELEN, args[1]));
		    obj.run(args[2], args[4]);
		}else{
		    System.err.print("USAGE:\t");
		    Grasper.printSVUsage();
		}
	    }else if(command.equals("sort")){
		if(args.length == 3){
		    SAMSort.main(Arrays.copyOfRange(args,1,3));
		}else{
		    System.err.print("USAGE:\t");
		    Grasper.printSortUsage();
		}
	    }else{
		System.err.println("Invalid command '" + command + "'");
		Grasper.printVersion();
	    }
	}else{
	    Grasper.printVersion();
	}
    }

    private static void printProgramList(){
	printDepthUsage();
	printSortUsage();
	printSVUsage();
    }

    private static void printDepthUsage(){
	System.err.println("java -jar grasper.jar depth <threadFile> <sorted sam> <config_file>\n");
    }
    private static void printSortUsage(){
	System.err.println("java -jar grasper.jar sort <bwa paired SAM file> <config_file>\n");
    }
    private static void printSVUsage(){
	System.err.println("java -jar grasper.jar sv <threadFile> <midsorted SAM from 'sort' program> <config_file> <depth file from 'depth' program>\n");
    }

    private static void printVersion(){
	System.err.println("\nProgram:\tGRAPSER (Genome Rearrangement Analysis using Short Paired-End Reads)");
	System.err.println("Version:\t0.1");
	System.err.println("Contact:\tHeewook Lee <heewlee@indiana.edu>");
	System.err.println("\nUsage:\t\tjava -jar grasper.jar <command> [options]\n");
	System.err.println("Command:\tdepth\tLoad read count (depth) information on edges from reads (serialization)");
	System.err.println("        \tsort\tRemoves obvious concordant read pairs and sort SAM files based on midpoint");
	System.err.println("        \tsv\tRemoves concordant read pairs, perform breakpoint clustering, and assign SV events\n");
    }
    
    private static boolean checkProgram(String command){
	if(command.equals("depth"))
	    return true;
	else if(command.equals("sv"))
	    return true;
	else if(command.equals("sort"))
	    return true;
	//else if(command.equals("full"))//legacy mode
	//    return true;
	else
	    return false;
    }
}
