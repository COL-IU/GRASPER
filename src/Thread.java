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

// Thread has genomeLength and intervalSequence as its classfield.
// intervalSequence is a ordered array of intervals.
// 
public class Thread{
    
    private int genomeLength;
    private ArrayList<Interval> intervalSequence;//sequence of Intervals representing the thread.
    //private Hashtable<String, ArrayList<Integer>> edgeId2IntervalListHash;

    
    public static void main(String[] args){
	//new Thread(4639675, args[0]).traverseDemo(); //ecoli.fna
	//new Thread(3284156, args[0]).traverseDemo(); //DeinoAll.fna
	new Thread(Constants.GENOMELEN, args[0]).traverseDemo(); //DeinoAll.fna
    }

    public Thread(int gLen, String threadFile){
	this.genomeLength = gLen;
	this.intervalSequence = new ArrayList<Interval>();
	new ThreadLoader().loadThread(threadFile, this);
    }
    
    public ArrayList<DirectionalRange> findDepthDepletedRegions(){
	ArrayList<DirectionalRange> deletedRanges = new ArrayList<DirectionalRange>();
	Interval startInterval = null;
	int startDepthArrIndex = -1;
	Interval prevInterval = null;
	int prevDepthArrIndex = -1;
	int count = 0;
	for(int k=0; k<this.intervalSequence.size(); k++){
	    Interval curInterval = this.intervalSequence.get(k);
	    boolean curIntDir = curInterval.isFwd();
	    Edge curEdge = curInterval.getEdge();
	    
	    if(!curInterval.getEdge().isRepetitive()){
		for(int i=0; i<curEdge.getDepthArrFwd().length;i++){
		    if(curEdge.getDepthAtIndex(i) < Constants.MIN_DEPTH){
			if(startInterval == null){
			    startInterval = curInterval;
			    startDepthArrIndex = i;
			    count++;
			}else{
			    count++;
			    prevInterval = curInterval;
			    prevDepthArrIndex = i;
			}
		    }else{//deletion ENDS and need to commit
			if(count >= 3){
			    deletedRanges.add(new DirectionalRange(startInterval.getLinPosFromStartDepthArrIndex(startDepthArrIndex)
								   , prevInterval.getLinPosFromEndDepthArrIndex(prevDepthArrIndex)
								   , true
								   , startInterval.getThreadIndex()
								   , prevInterval.getThreadIndex()));
			}
			startInterval = null;
			startDepthArrIndex = -1;
			prevInterval = null;
			prevDepthArrIndex = -1;
			count = 0;
		    }
		}
	    }
	}
	return deletedRanges;
    }


    public void serializeDepthArray(String samfile){
	System.err.println("*** Serializing Depth Array . . .");
	//BufferedWriter bw = null;
	FileOutputStream fos = null;
	ObjectOutputStream oos = null;
	//before serializing depth array, we wrap it around using ArrayList. then serialize the entire ArrayList
	ArrayList<int[]> list = new ArrayList<int[]>();
	for(int i=0; i<this.intervalSequence.size();i++){
	    Interval curInterval = this.intervalSequence.get(i);
	    if(curInterval.getEC() == 0){
		if(Constants.DEBUG3){
		    int[] edgeLackingValArr = curInterval.getEdge().isEdgeLackingCoverage();
		    System.err.println("[DEPTH][" + curInterval.getEdge().getId() + "]\t\t" 
				       + edgeLackingValArr[0] + ":" +edgeLackingValArr[1] + ":" + edgeLackingValArr[2]);
		}
		list.add(curInterval.getEdge().getDepthArrFwd());
		list.add(curInterval.getEdge().getDepthArrRev());
	    }
	}
	
	try{
	    //bw = new BufferedWriter(new FileWriter(samfile+".depth"));
	    fos = new FileOutputStream(samfile+".depth");
	    oos = new ObjectOutputStream(fos);
	    oos.writeObject(list);
	    oos.close();
	    fos.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
    }

    public void deserializeDepthArray(String data){
	
	ArrayList<int[]> list = new ArrayList<int[]>();
	try{
	    FileInputStream fis = new FileInputStream(data);
	    ObjectInputStream ois = new ObjectInputStream(fis);
	    list = (ArrayList<int[]>) ois.readObject();
	    ois.close();
	    fis.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	    return;
	}catch(ClassNotFoundException c){
	    System.err.println("De-Serializer: Class not found");
	    c.printStackTrace();
	    return;
	}

	int count = 0;
	for(int i=0; i<this.intervalSequence.size();i++){
	    Interval curInterval = this.intervalSequence.get(i);
	    if(curInterval.getEC() == 0){
		curInterval.getEdge().setDepthArrFwd(list.get(count));
		count++;
		curInterval.getEdge().setDepthArrRev(list.get(count));
		count++;
		if(Constants.DEBUG3){
		    int[] edgeLackingValArr = curInterval.getEdge().isEdgeLackingCoverage();
		    System.err.println("[DEPTH][" + curInterval.getEdge().getId() + "]\t\t" 
				       + edgeLackingValArr[0] + ":" +edgeLackingValArr[1] + ":" + edgeLackingValArr[2]);
		}
	    }
	}
	
	System.err.println("====== Done De-Serializing");
    }

    public int[] doesNOverlapOnRange(int n, int s, int e){
	
	GraphPos maxGP = this.getGraphPosForPositionX(n);
	int[] linPoss = maxGP.getEdge().getLinearPositions(maxGP.getGraphPos());
	
	ArrayList<Integer> overlapPositions = new ArrayList<Integer>();
	for(int i=0; i<linPoss.length; i++){
	    if( linPoss[i] >= s 
		&& linPoss[i] <= e )
		overlapPositions.add(new Integer(linPoss[i] - s + 1));
	}

	int[] poss = new int[overlapPositions.size()];
	for(int i=0;i<overlapPositions.size(); i++){
	    poss[i] = overlapPositions.get(i).intValue();
	}
	return poss;
    }



    /*
    public double ratioOfInsufficientDepthBins(Edge e1, int gp1, int ec1, Edge e2, int gp2, int ec2){
	Edge edge1 = null;
	int graphPos1 = gp1;
	int edgeCounter1 = ec1;
	Edge edge2 = null;
	int graphPos2 = gp2;
	int edgeCounter2 = ec2;
	
	if(e1.nthInterval(ec1) >)
	    ;
	
	    }*/

    /*
     * returns :
     * 1 : if the path lacks coverage
     * 0 : if the path has coverages and edge has to be non-repetative
     * -1 : if the path has coverages and edge has to be repetative.
     */
    public int isPathLackingCoverage(IntervalPos ip1, IntervalPos ip2){
	return this.isPathLackingCoverage(ip1.getInterval(), ip1.getGraphPos(), ip2.getInterval(), ip2.getGraphPos());
    }
    
    
    /*
     * returns :
     * 1 : if the path lacks coverage
     * 0 : if the path has coverages and edge has to be non-repetative
     * -1 : if the path has coverages and edge has to be repetative.
     */
    public int isPathLackingCoverage(Interval i1, int gp1, Interval i2, int gp2){
	
	Interval itv1 = i1;
	int graphPos1 = gp1;
	Interval itv2 = i2;
	int graphPos2 = gp2;
	
	if(i1.getThreadIndex() > i2.getThreadIndex()){
	    itv1 = i2;
	    graphPos1 = gp2;
	    itv2 = i1;
	    graphPos2 = gp1;
	}else if(i1.getThreadIndex() == i2.getThreadIndex()){
	    if(gp1 > gp2){
		graphPos1 = gp2;
		graphPos2 = gp1;
	    }
	}

	int tIndex1 = itv1.getThreadIndex();
	int tIndex2 = itv2.getThreadIndex();
	int diff = tIndex2 - tIndex1;

	int pathLength = 0; //this included the entire length of the path
	int totalLength = 0; //this is the total lengths of segments we are actually counting. totalLegnth <= pathLength
	int lackingLength = 0; //this is the lengths of segments that are deemed lacking coverages. lackingLength <= pathLength
	
	//if on same interval
	if(diff == 0){
	    int lackingVal = itv1.getEdge().isPartialEdgeLackingCoverage(graphPos1, graphPos2);
	    if(lackingVal == 1){//1 means --> edge is lacking coverages, regardless of whether Edge is repetative or not
		return 1;//return true;
	    }else if(lackingVal == -1){ //edge has coverage but it's repetative.
		
		return -1;
	    }else
		return 0;
	    //lackingLength = (itv1.getEdge().isPartialEdgeLackingCoverage(graphPos1, graphPos2) ? (graphPos2-graphPos1+1) : 0);
	    //totalLength = graphPos2-graphPos1+1;
	}else if(diff > 0){
	    //first interval (partial : from gp to end of interval)
	    //lackingValArr [0]: lackingVal,  [1]: insuffDepthbinsCount,  [2]: totalBinsCount
	    int[] lackingValArrFirst = itv1.getEdge().isPartialEdgeLackingCoverage(graphPos1, itv1.isFwd(), false);
	    /*pathLength = itv1.getEdge().computeSegmentLength(graphPos1, itv1.isFwd(), false);
	    if(lackingVal == 1){
		lackingLength = pathLength;
		totalLength = pathLength;
	    }else if(lackingVal == 0)
		totalLength = pathLength;
	    */

	    //last interval (partial : from beginning of interval to gp)
	    int[] lackingValArrLast = itv2.getEdge().isPartialEdgeLackingCoverage(graphPos2, itv2.isFwd(), true);
	    //int tmpLen = itv2.getEdge().computeSegmentLength(graphPos2, itv2.isFwd(), true);
	    /*pathLength += tmpLen;
	    if(lackingVal == 1){
		lackingLength += tmpLen;
		totalLength += tmpLen;
	    }else if(lackingVal == 0)
		totalLength += tmpLen;
	    */
	    //intervals in between
	    int[] lackingValArrMiddle = new int[3];
	    for(int i=(tIndex1+1); i<tIndex2; i++){
		Interval itv = this.intervalSequence.get(i);
		int[]lackingValArrTmp = itv.getEdge().isEdgeLackingCoverage();
		//pathLength += itv.getEdge().length();
		if(lackingValArrTmp[0] != -1){
		    lackingValArrMiddle[1] += lackingValArrTmp[1];
		    lackingValArrMiddle[2] += lackingValArrTmp[2];
		}
	    }
	    int[] lackingValTotal = new int[2]; // 0: insufficientBinsCount, 1: totalBinsCount
	    lackingValTotal[0] = (lackingValArrFirst[0] == -1 ? 0 : lackingValArrFirst[1] ) 
		+ (lackingValArrLast[0] == -1 ? 0 : lackingValArrLast[1] ) 
		+ lackingValArrMiddle[1];
	    lackingValTotal[1] = (lackingValArrFirst[0] == -1 ? 0 : lackingValArrFirst[2] ) 
		+ (lackingValArrLast[0] == -1 ? 0 : lackingValArrLast[2] ) 
		+ lackingValArrMiddle[2];
	    if(lackingValTotal[1] == 0)//repetative or too short
		return -1;
	    else{
		double lackingBinsRatio = (lackingValTotal[0]*1.0d) / (lackingValTotal[1]*1.0d);
		if(lackingBinsRatio > Constants.NO_COVERAGE_RATIO)
		    return 1;
		else
		    return 0;
	    }
	    
	}else
	    return -5;//returns something else. it never reaches here.

    }
    
    		


    public ArrayList<Interval> getIntervalSequence(){
	return this.intervalSequence;
    }

    public void traverseDemo(){
	int lenSum = 0;
	System.err.println("=======================================");
	System.err.println("= genomeLength : " + this.genomeLength);
	System.err.println("= numIntervals : " + this.intervalSequence.size());
	System.err.println("=======================================");
	System.err.println();
	//int tmp = 1;
	for(int i=0;i<intervalSequence.size();i++){
	    System.err.println("@" + i + "\t-->\t" + intervalSequence.get(i).toString());
	    //if(tmp-1 != intervalSequence.get(i).getStart())
	    //		System.out.println("---------> CHECK!");
	    //tmp = intervalSequence.get(i).getStart() + intervalSequence.get(i).length();
	    lenSum += intervalSequence.get(i).length();
	}
	System.err.println();
	System.err.println("=======================================");
	System.err.println("= lenSum : " + lenSum );
	//this.searchDemo(1000000);
    }

    public void searchDemo(int numSearches){
	int[] posArr = new int[numSearches];
	Random generator = new Random();
	for(int i=0;i<posArr.length;i++)
	    posArr[i] = generator.nextInt(Constants.GENOMELEN) + 1;
	    //posArr[i] = generator.nextInt(4639675) + 1;
	
	for(int i=0; i<posArr.length;i++)
	    System.out.println("Query(x): " + posArr[i] + "\t==>\t"+this.getIntervalForPositionX(posArr[i]));
	
    }

    public int length(){
	return this.genomeLength;
    }
    
    public int getNumIntervals(){
	return this.intervalSequence.size();
    }

    public void add(Interval i){
	this.intervalSequence.add(i);
    }
    
    public Interval getNthInterval(int n){
	return this.intervalSequence.get(n);
    }

    public Edge getEdgeFromNthInterval(int n){
	return this.getNthInterval(n).getEdge();
    }
    
    /* binary search of the thread to find the interval for position x. x is linear position on the genome*/
    /*
    public Interval getIntervalForPositionX(int x, int n){
	int maxIndex = this.intervalSequence.size()-1;
	int minIndex = 0;
	int curIndex = maxIndex/2;
	boolean contin = true;
	int compVal = 0;
	while(contin){
	    //System.err.print(curIndex + "[" + minIndex + "," + maxIndex +"]");
	    compVal = this.intervalSequence.get(curIndex).compare(x);
	    if(compVal < 0){//if position x is before the interval
		maxIndex = curIndex - 1;
		curIndex = curIndex - (curIndex- minIndex+1)/2;
	    }else if(compVal > 0){
		minIndex = curIndex + 1;
		curIndex = curIndex + (maxIndex - curIndex+1)/2;
	    }else
		break;
	}
	return this.intervalSequence.get(curIndex);
	}*/
    
    /* binary search of the thread to find the interval for position x. x is linear position on the genome*/
    /* returns null if can't find */
    public Interval getIntervalForPositionX(int x){
	int maxIndex = this.intervalSequence.size()-1;
	int minIndex = 0;
	while(maxIndex >= minIndex){
	    int curIndex = (minIndex + maxIndex)/2; 
	    int compVal = this.intervalSequence.get(curIndex).compare(x);
	    if(compVal < 0)
		maxIndex = curIndex - 1;
	    else if(compVal > 0)
		minIndex = curIndex + 1;
	    else
		return this.intervalSequence.get(curIndex);
	}
	
	return null;
    }

    public GraphPos getGraphPosForPositionX(int x){
	//Interval intv = this.getIntervalForPositionX(x);
	//int graphPos = intv.linPos2GraphPos(x);
	//return new GraphPos(intv.getEdge(), graphPos);
	return new GraphPos(this, x);
    }

    public IntervalPos getIntervalPosForPositionX(int x){
	//Interval intv = this.getIntervalForPositionX(x);
	//int intPos = intv.linPos2IntervalPos(x);
	//return IntervalPos(intv, intPos);
	return new IntervalPos(this, x);
    }
    
    public int[] getIntervalsForRange(int s, int e){
	
	int startTI = this.getIntervalForPositionX(s).getThreadIndex();
	int endTI = this.getIntervalForPositionX(e).getThreadIndex();
	int[] tis = new int[endTI - startTI + 1];
	
	for(int i=0; i<tis.length; i++)
	    tis[i] = startTI + i;
	
	return  tis;
    }

    
}

class IntervalPos{

    private Interval intv;
    private int gp;

    public Interval getInterval(){
	return this.intv;
    }
    
    public int getGraphPos(){
	return this.gp;
    }

    public int getThreadIndex(){
	return this.intv.getThreadIndex();
    }
        
    public IntervalPos(Interval i, int p){
	this.intv = i;
	this.gp = p;
    }

    public IntervalPos(Thread t, int linPos){
	this.intv = t.getIntervalForPositionX(linPos);
	this.gp = this.intv.linPos2GraphPos(linPos);
    }
    
}

class GraphPos{
    private Edge edge;
    private int gp;
    
    public Edge getEdge(){
	return this.edge;
    }
    
    public int getGraphPos(){
	return this.gp;
    }
    
    public GraphPos(Edge e, int p){
	this.edge = e;
	this.gp = p;
    }
    
    public GraphPos(Thread t, int linPos){
	Interval tmp = t.getIntervalForPositionX(linPos);
	this.edge = tmp.getEdge();
	this.gp = tmp.linPos2GraphPos(linPos);
    }
    
}

class ThreadLoader{
    
    private HashMap<String, Interval> intervalHash;
    
    public ThreadLoader(){
	intervalHash = new HashMap<String, Interval>();
    }
    
    public Thread loadThread(String threadFile, Thread t){
	BufferedReader br = null;
	try{
	    br = new BufferedReader(new FileReader(threadFile));
	    String curline = "";
	    //int ec = 0;
	    int curThreadIndex = 0;
	    while( (curline=br.readLine()) != null ){
		if(!curline.startsWith("#")){
		    if(Constants.DEBUG)
			System.err.println("processing:\t" + curline);
		    t.add(processInterval(curline, t.length(), curThreadIndex));
		    curThreadIndex++;
		}
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	return t;
    }
    
    /* threadFile interval Line format:
       field 0 : edgeID1
       field 1 : edgeID2
       field 2 : multiplicity
       field 3 : edge length
       
       field 4i :   direction --> 0:forwad, direction 1:reverse
       field 4i+1 : start  
       field 4i+2 : length
    */
    private Interval processInterval(String intervalStr, int genomeLength, int threadIndex){
	String[] tokens = intervalStr.split("\\t");
	return this.processInterval(Integer.parseInt(tokens[0]), 
				    Integer.parseInt(tokens[1]), 
				    Integer.parseInt(tokens[2]), 
				    Integer.parseInt(tokens[3]), 
				    tokens,
				    genomeLength, 
				    threadIndex);
    }
    
    
    private Interval processInterval(int ei1, int ei2, int m, int l,  String[] tokens, int genomeLength, int threadIndex){
	String curKey = this.getKeyFromEdgeIDs(ei1,ei2);
	if(Constants.DEBUG)
	    System.out.println("Loading the edge : [" + curKey+"]" );
	Interval tmpInterval = null;

	/*have we already visited the edge meaning, already in the hash?*/
	
	//if NOT, we need to create the new edge and update.
	if( (tmpInterval = this.intervalHash.remove(curKey)) == null){
	    
	    Edge tmpEdge = new Edge(curKey , m , tokens , genomeLength);
	    tmpInterval = tmpEdge.nthInterval(0);//since we are newly adding this edge, this must be the 0-th interval.
	    
	}else{
	    //System.out.println(tmpInterval);
	    tmpInterval = tmpInterval.getNextInterval(); //if already in the hash, need to get a new instance of interval with updated values.
	}

	tmpInterval.setThreadIndex(threadIndex); //setting thread index value to interval.

	/* need to put it into the hash only if the current occurence of the edge (current interval) is 
	   not the last occurence */
	if(!tmpInterval.isLastInterval())
	    this.intervalHash.put(curKey, tmpInterval);
	
	return tmpInterval;
    }
    
    //key is always -->  smaller number : larger number
    private String getKeyFromEdgeIDs(int ei1, int ei2){
	if(ei1 < ei2)
	    return ei1 + ":" + ei2;
	else
	    return ei2 + ":" + ei1;
    }
    
}
