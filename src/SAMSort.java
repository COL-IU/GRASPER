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
import java.util.*;

public class SAMSort{

    public static void main(String[] args) throws Exception{

	if(args.length != 2){
	    System.err.println("USAGE: java -jar SAMSort.jar <BWA PAIRED SAM file> <configFile>");
	}else{
	    System.err.println("Running: java -jar SAMSort.jar " + args[0] + " " + args[1]);
	    //args[1]: conf file
	    Constants.loadConstants(args[1], false);
	    
	    StringBuffer headers = new StringBuffer();
	    
	    ArrayList<SAM> list = new ArrayList<SAM>();
	    
	    BufferedReader br = new BufferedReader(new FileReader(args[0]));
	    int i = 0;
	    int count = 0;
	    String curline = null;
	    String curline2 = null;
	    
	    SAM s1 = null;
	    SAM s2 = null;
	    
	    boolean done = false;
	    
	    //int cordantCount = 0;
	    //int discordantCount = 0;
	    
	    int totalPairCount = 0;
	    int discordantPairCount = 0;
	    int singleEndMappedPairCount = 0;
	    int unmappedPairCount = 0;

	    BufferedWriter unmapped_bw = null;
	    BufferedWriter se_bw = null;

	    while( !done ){
		curline = br.readLine();
		if(curline == null)
		    break;
		else if(!curline.startsWith("@")){
		    s1 = new SAM(curline);
		    //skip secondary. Load until primary
		    while(s1.isSecondary()){
			curline = br.readLine();
			if(curline == null){
			    done = true;
			    break;
			}
			s1 = new SAM(curline);
		    }
		    
		    if(done)
			break;
		    
		    curline2 = br.readLine();
		    s2 = new SAM(curline2);
		    //skip secondary. Load until primary
		    while(s2.isSecondary()){
			curline = br.readLine();
			if(curline == null){
			    done = true;
			    break;
			}
			s2 = new SAM(curline);
		    }
		    
		    if(done)
			break;
		    
		    totalPairCount++;
		    

		    //only collect discordant
		    SAMPair tmppair = new SAMPair(s1, s2);
		    if(tmppair.isDiscordant()){
			discordantPairCount++;
			list.add(s1);
			list.add(s2);
			s1 = null;
			s2 = null;
			if(list.size() == 200000){
			    Collections.sort(list);
			    BufferedWriter bw = new BufferedWriter(new FileWriter(args[0] + ".part_" + i));
			    for( SAM s: list){
				bw.write(s.getSamline() + "\n");
			    }
			    
			    bw.close();
			    bw = null;
			    i++;
			    list = new ArrayList<SAM>();
			}
		    }else if(!tmppair.isBothUnmapped()){//single-end mapped
			if(se_bw == null){
			    se_bw = new BufferedWriter(new FileWriter(args[0] + ".singleEndMapped"));
			    se_bw.write(headers.toString());
			}
			singleEndMappedPairCount++;
			se_bw.write(s1.getSamline() + "\n");
			se_bw.write(s2.getSamline() + "\n");
			s1 = null;
			s2 = null;
		    }else{ // both unmapped
			if(unmapped_bw == null){
			    unmapped_bw = new BufferedWriter(new FileWriter(args[0] + ".unmapped"));
			    unmapped_bw.write(headers.toString());
			}
			unmappedPairCount++;
			unmapped_bw.write(s1.getSamline() + "\n");
			unmapped_bw.write(s2.getSamline() + "\n");
			s1 = null;
			s2 = null;
		    }
		    
		}else{
		    headers.append(curline + "\n");
		}
		
	    }
	    br.close();
	    br = null;
	    
	    if(se_bw != null)
		se_bw.close();
	    if(unmapped_bw != null)
		unmapped_bw.close();

	    if(list.size() > 0){
		Collections.sort(list);
		BufferedWriter bw = new BufferedWriter(new FileWriter(args[0] + ".part_" + i));
		//System.err.println(count + "\t" + set.size());
		for(SAM s: list){
		    bw.write(s.getSamline() + "\n");
		}
		bw.close();
		bw = null;
		i++;
		list = null;
	    }
	    
	    Map<SAM, Queue<SAM>> map = new TreeMap<SAM, Queue<SAM>>();
	    //Set<SAM> set = new TreeSet<SAM>();
	    //System.err.println("openning " + i + " files");
	    BufferedReader[] brArr = new BufferedReader[i];
	    for(int j=0; j<i; j++){
		brArr[j] = new BufferedReader(new FileReader(args[0] + ".part_" + j));
		String tmpline = brArr[j].readLine();
		SAM tmp = new SAM(tmpline, j);
		if(!map.containsKey(tmp))
		    map.put(tmp, new LinkedList<SAM>());
		else{
		    //System.err.println("here");
		    map.get(tmp).add(tmp);//brArr[j].pushBack(tmpline);
		}
		//System.err.println(map.get(tmp).size());
	    }
	    //System.err.println(map.size());
	    
	    BufferedWriter bw = new BufferedWriter(new FileWriter(args[0]+".discordant.midsorted"));
	    bw.write(headers.toString());
	    SAM s = null;
	    
	    //int counter = 0;
	    Queue<SAM> q = null;
	    while(!map.isEmpty()){
		s = map.keySet().iterator().next();
		q = map.get(s);
		//map.remove(s);
		bw.write(s.getSamline());
		bw.write("\n");
		if(q.size() == 0){
		    map.remove(s);
		    int tmpIndex = s.getIndex();
		    String tmpline = brArr[tmpIndex].readLine();
		    if(tmpline != null){
			s = new SAM(tmpline, tmpIndex);
			if(!map.containsKey(s)){
			    //System.err.print("put b4: " + map.size());
			    map.put(s, new LinkedList<SAM>());
			    //System.err.println("\ta4: " + map.size());
			}else{
			    //System.err.print("here2: ");
			    map.get(s).add(s);
			    //System.err.println(map.size());
			}
		    }//else{
		    //System.err.println("reader[" + tmpIndex + "] : NULL");
		    //}
		}else{
		    //System.err.print("B4 Q size: " + q.size() + " map size: " + map.size() );
		    map.remove(s);
		    map.put(q.remove(), q);
		    //System.err.println("\tA4 Q size: " + q.size() + " map size: " + map.size() );
		    int tmpIndex = s.getIndex();
		    String tmpline = brArr[tmpIndex].readLine();
		    if(tmpline != null){
			s = new SAM(tmpline, tmpIndex);
			if(!map.containsKey(s)){
			    //System.err.print("put b4: " + map.size());
			    map.put(s, new LinkedList<SAM>());
			    //System.err.println("\ta4: " + map.size());
			}else{
			    //System.err.print("here2: ");
			    map.get(s).add(s);
			    //System.err.println(map.size());
			}
		    }//else{
		    //  System.err.println("reader[" + tmpIndex + "] : NULL");
		    //}
		    
		    
		}
		//counter++;
		
	    }
	    //System.err.println(counter);
	    bw.close();
	    
	    for(int j=0;j<brArr.length;j++){
		brArr[j].close();
		new File(args[0] + ".part_" + j).delete();
	    }
	    
	    //	int totalPairCount = 0;
	    //int discordantPairCount = 0;
	    System.out.println("------ DONE WITH REMOVAL OF CONCORDANT & UNMAPPED PAIRS ---------");
	    System.out.println("------ Processed a total of\t" + totalPairCount + " read-pairs" );
	    System.out.println("------ Retained a total of\t" + discordantPairCount + " read-pairs writtend out to " + args[0] + ".discordant.midsorted" );
	    System.out.println("------ Single-end mapped :\t" + singleEndMappedPairCount + " read-pairs written out to " + args[0] + ".singleEndMapped" );
	    System.out.println("------ Unmapped (both ends) : \t" + unmappedPairCount + " read-pairs written out to " + args[0] + ".unmapped" );
	}
    }
    
}

class BufferedReaderPB{
    
    private String lastLine; 
    private BufferedReader reader;
    
    public BufferedReaderPB(String file) throws Exception{
	this.lastLine = null;
	this.reader = new BufferedReader(new FileReader(file));
    }

    void close() throws Exception{
	this.reader.close();
    }

    String readLine() throws Exception{
	if(lastLine != null){
	    String tmp = lastLine;
	    lastLine = null;
	    return tmp;
	}else
	    return this.reader.readLine();
    }
    
    void pushBack(String line){
	this.lastLine = line;
    }
}

