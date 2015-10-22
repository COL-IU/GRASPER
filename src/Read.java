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
import java.util.regex.*;

/////////////////////////////////////////////////////////////////////////////////////////////
//    private Edge mappedEdge;
//    private int graphPos; //distance from 5' side of the forward traversal interval. Must flip for reversal. //UPDATED so that this is the midpoint of this read.
//    private boolean fwd; // NEED TO UPDATE THIS DIRECTION so that it's with the respect to the first interval of the mappedEdge
//    private Read pairingRead;
//    private int length;
//    private int linPos;
//    private String readName;
//    private String genome;
//    private boolean unique;
/////////////////////////////////////////////////////////////////////////////////////////////
public class Read implements Comparable<Read>{
    
    public Read(){
	//NEED TO WRITE sameline parser.
	this.mappedEdge = null;
	this.graphPos = -1;/* samline position : LINEAR GENOMIC POSITION! --> to graphPOS */
	this.fwd = false;  /* samline direction : Linear Reference direction */ // this is later UPDATED respect to the first interval
	this.pairingRead = null;
	this.length = 0; /* samline readlength */
	this.needToDelete = false;
	this.firstRead = true;
	
    }

    //initial constructor
    // this constructor only sets up the read information
    // This must be updated before the object is actually loaded onto the graph.
    public Read(String[] tokens){
	this();
	//this.linPos = Integer.parseInt(tokens[3]);
	this.length = this.getEffectiveLengthFromCIGAR(tokens[5]); //tokens[9].length();
	this.readName = tokens[0];
	this.genome = tokens[2];
	this.linPos = Integer.parseInt(tokens[3]); // left-most position in SAM
	this.fwd = getDirection(Short.parseShort(tokens[1]));
	this.updateLinPos();
	this.uniq = isUniqueMapping(tokens);
    }

    public boolean isReadRepetitiveOnUniqueEdge(){
	if(this.uniq)
	    return false;
	else{
	    if(this.mappedEdge.isUniqueEdge()){
		if(Constants.DEBUG2)
		    System.err.println("[Read.isReadRepetitiveOnUniqueEge] REMOVING THIS MAPPING BECAUSE OF MULTIPLE MAPPING on unique EDGE. :\t" + this.readName);

		return true;
	    }
	    return false;
	}
    }

    //we only consider M and I for the effective length of READ.
    public int getEffectiveLengthFromCIGAR(String cigar){
	Pattern pattern = Pattern.compile("([1-9]\\d*)([MIDNSHPX=])");
	Matcher matcher = pattern.matcher(cigar);
	int sumMI = 0;
	while(matcher.find()){
	    //int sumMI = 0;
	    int size = Integer.parseInt(matcher.group(1));
	    char type = matcher.group(2).charAt(0);
	    if(type == 'M' || type == 'I')
		sumMI += size;
	}
	return sumMI;
    }

    private boolean isUniqueMapping(String[] tokens){
	try{
	    int alignscore = Integer.parseInt(tokens[13].substring(tokens[13].lastIndexOf(":")+1));
	    int suboptimalscore = Integer.parseInt(tokens[14].substring(tokens[13].lastIndexOf(":")+1));
	    if(alignscore > suboptimalscore)
		return true;
	}catch(Exception e){
	    e.printStackTrace();
	    for(int i=0;i<tokens.length;i++)
		System.err.print("\t" + tokens[i]);
	    System.err.println();
	    System.exit(0);
	}
	return false;
    }



    //updates the linpos to midpost rather than 5' left position in linear coordinate
    private void updateLinPos(){
	this.linPos = this.linPos + this.getDistanceToEnd(); //
    }

    public String toString(){
	//return readName + "\t" + linPos + "\t" + fwd + "\tEdgeID:" + mappedEdge.getId() +"\t"+(pairingRead != null ? pairingRead.getReadName() : "null");
	return readName + "\t" + this.metaString() + "<||>" + (pairingRead != null ? pairingRead.metaString() : "null");
    }

    public String toStringAtEc(int ec){
	return (this.firstRead ? "(1)" : "(2)") + readName + "\t" + this.metaStringAtEc(ec) + "<||>" + (pairingRead != null ? pairingRead.metaString() : "null");
    }

    public String metaString(){
	//return linPos + "(" + (fwd ? "+" : "-") + ")|" + mappedEdge.getId()+"|graphPos:"+this.graphPos;
	String tmp = null;
	try{
	    tmp = "SAMmidPos:" + linPos + "(" + (fwd ? "+" : "-") + ")|" + mappedEdge.getId()+"|graphPos:"+this.graphPos;
	}catch(NullPointerException npe){
	    System.err.println("ERROR CAUGHT:\t" + this.readName);
	    npe.printStackTrace();
	    return null;//System.exit(0);
	}
	return tmp;
    }

    public String metaStringAtEc(int ec){
	String tmp = null;
	try{
	    //tmp = linPos + "(" + (fwd ? "+" : "-") + ")|" + mappedEdge.getId()+"|graphPos:"+this.graphPos;
	    tmp = "curLinPos:" + this.getLinPosWithEdgeCounterAt(ec) + "("+ (this.isFwdAtEc(ec) ? "+" : "-") +")|" + mappedEdge.getId()+"|graphPos:"+this.graphPos;
	}catch(NullPointerException npe){
	    System.err.println("ERROR CAUGHT:\t" + this.readName);
	    npe.printStackTrace();
	    return null;//System.exit(0);
	}
	return tmp;
    }

    
    private static boolean getDirection(short flag){
	short reverse = 0x0010;
	if((flag & reverse) == 0)
	    return true;
	return false;
    }
    

    //no longer used : merged into setGraphPosAndDirection(Interval)
    private void updateDirection(Interval intv){
	if(!intv.isFwd())
	    this.fwd = !this.fwd;
    }

    public void update(Interval intv){
	this.mappedEdge = intv.getEdge();
	//this.linPos = 
	this.setGraphPosAndDirection(intv);
	//this.updateDirection(intv); /*merged into setGraphPosAndDirection(Interval)*/
    }
    
    public void updateDepth(Interval intv){
	this.mappedEdge = intv.getEdge();
	//this.linPos = 
	this.setGraphPosAndDirectionDepth(intv);
	//this.updateDirection(intv); /*merged into setGraphPosAndDirection(Interval)*/
    }

    

    /* ONLY INVOKE THIS IF TWO READS ARE BOTH MAPPED to same EDGE*/
    /*
    private int distanceToWithinSameEdge(Read other){
	if(this.graphPos > other.getGraphPos()){
	    if( (this.graphPos + this.getDistanceToEnd() ) > ( other.getGraphPos() + other.getDistanceToEnd() ) )
		return this.graphPos + this.getDistanceToEnd() - other.getGraphPos() + other.getDistanceToEnd();
	    else
		return other.getLength();
	}else{
	    if( ( this.graphPos + this.getDistanceToEnd() ) < (other.getGraphPos() + other.getDistanceToEnd() ) )
		return other.getGraphPos() + other.getDistanceToEnd() - this.graphPos + this.getDistanceToEnd();
	    else
		return this.length;
	}
    }
    */
    
    //new version: we only consider distance between midpoints (graphPos)
    //assuming reads are approximately same length
    private int distanceToWithinSameEdge(Read other){
	if(this.graphPos > other.getGraphPos())
	    return this.graphPos - other.getGraphPos() + 1;
	else
	    return other.getGraphPos() - this.graphPos + 1;
    }

    public int getDistanceToEnd(){
	return this.length/2;
    }

    private boolean isFromSameEdge(Read other){
	if(this.mappedEdge.equals(other.getMappedEdge()))
	    return true;
	return false;
    }

    public int distanceTo(Read other, int maxDistance){
	//if from same edges
	if(this.isFromSameEdge(other))
	    return this.distanceToWithinSameEdge(other);
	else
	    return this.mappedEdge.getDistanceFromGraphPosTo(this.graphPos, other, maxDistance);
    }

    /*
    public int distanceTo(Read other){
	if(this.mappedEdge.equals(other.getMappedEdge()))
	    return Math.abs(this.graphPos - other.graphPos) + length;
	else
	    return this.mappedEdge.getDistanceFromGraphPosTo(this.graphPos, other);
    }
    */
    
    
    //assuming the read length is same.
    //this returns distance between currentLinearPosition and a read(other)
    //if the read(other) is on the same edge, simply calculate the distance.
    //else --> find two closest intervals
    public int distanceToOtherReadFromCurrentPosition(int ec, Read other){
	if(this.isFromSameEdge(other))
	    return this.distanceToWithinSameEdge(other);
	else
	    return other.getMappedEdge().getDistanceToClosestInterval(other.getGraphPos(), this.getLinPos() );
    }


    // NEED TO UPDATE TO HANDLE READS ON OTHER EDGES
    // NOT USED CURRENTLY.
    public int compareTo(Read other){
	if(other == null)
	    return -1;
	if(this.graphPos < other.getGraphPos())
	    return -1;
	else if(this.graphPos == other.getGraphPos())
	    return 0;
	else
	    return 1;
    }
    
    public int compareToReadOnSameEdge(Read other){
	if(other == null)
	    return -1;
	
	if(this.graphPos < other.getGraphPos())
	    return -1;
	else if(this.graphPos == other.getGraphPos())
	    return 0;
	else
	    return 1;
    }

    //compare this read to other read that are on same edge at ec-th interval.
    public boolean isBeforeOnSameEdgeAtEc(Read other, int ec){
	if(other == null){
	    if(Constants.DEBUG)
		System.err.println("WTF!!");
	}
	//System.err.println("Mapped Edge numIntervals:\t" + this.mappedEdge.getMultiplicity() + "\t,\t passedEC:\t" + ec );
	if(this.mappedEdge.nthInterval(ec).isFwd()){
	    if(this.compareToReadOnSameEdge(other) <= 0)
		return true;
	    else
		return false;
	}else{
	    if(this.compareToReadOnSameEdge(other) >= 0)
		return true;
	    else
		return false;
	}
	    
    }


    
    /////---------------CHECK--------------///////
    public int getLinPosWithEdgeCounterAt(int ec){
	return this.mappedEdge.nthInterval(ec).getLinearPosition(this.graphPos);
    }

    /* returns read's directon on ec-th interval. The direction is genome direction */
    public boolean isFwdAtEc(int ec){
	//return ( this.mappedEdge.nthInterval(ec).isFwd() && this.isFwd() );
	return ( this.mappedEdge.nthInterval(ec).isFwd() == this.isFwd() );
    }

    
    //---------------------------SETTERS----------------------------//

    /////---------------CHECK--------------///////
    public void setGraphPos(Interval segment){
	this.mappedEdge = segment.getEdge();
	//each interval belonging to the same edge can have different(similar) lengths. So here we normalize the length so that it fits within the edge length(length of the first segment)
	this.graphPos = (int) ( ((this.linPos - segment.getStart() + 1)*1.0d)
				/ (1.0d*this.mappedEdge.nthInterval(segment.getEC()).length())
				* (this.mappedEdge.length()*1.0d )
				);
	//this.graphPos - segment.getStart() + 1; //subtract start location of the segment linear genomic position 
	/* sets the direction up correctly */
	if(!this.mappedEdge.getIntervals()[segment.getEC()].isFwd()){
	    this.graphPos = this.mappedEdge.length() - this.graphPos;
	    //this.fwd = !(this.fwd); //done in updateDirection
	}
    }
    
    private void setGraphPosAndDirectionDepth(Interval segment){
	if(Constants.DEBUG4)
	    System.err.println("\t[setGraphPosAndDirection]b4:\tlinPos[" + this.linPos + "]\t-->\tgraphPos["+this.graphPos+"]" );
	this.graphPos = segment.linPos2GraphPos(this.linPos);
	updateDirection(segment);
	if(Constants.DEBUG4)
	    System.err.println("\t[setGraphPosAndDirection]after:\tlinPos[" + this.linPos + "]\t-->\tgraphPos["+this.graphPos+"]" );
	segment.getEdge().incrementDepthAt(this.graphPos, this.fwd); //this increments bin depth on EDGE.
    }

    private void setGraphPosAndDirection(Interval segment){
	if(Constants.DEBUG4)
	    System.err.println("\t[setGraphPosAndDirection]b4:\tlinPos[" + this.linPos + "]\t-->\tgraphPos["+this.graphPos+"]" );
	this.graphPos = segment.linPos2GraphPos(this.linPos);
	updateDirection(segment);
	if(Constants.DEBUG4)
	    System.err.println("\t[setGraphPosAndDirection]after:\tlinPos[" + this.linPos + "]\t-->\tgraphPos["+this.graphPos+"]" );
	//segment.getEdge().incrementDepthAt(this.graphPos, this.fwd); //this increments bin depth on EDGE.
    }

    //all this is moved to linPos2GraphPos(int linPos) in Interval.java
    /*
    private void setGraphPosAndDirection(Interval segment){
	this.mappedEdge = segment.getEdge();
	//each interval belonging to the same edge can have different(similar) lengths. So here we normalize the length so that it fits within the edge length(length of the first segment)
	this.graphPos = (int) ( ((this.linPos - segment.getStart() + 1)*1.0d)
				/ (1.0d * segment.length())
				* (this.mappedEdge.length()*1.0d) 
				);
	// for reverse segment, gp has to be flipped. No need to be offset by half-readLength as linPos used to compute gp is already midpoint
	// initially read direction is loaded according to SAM outputs (that is just respect to the reference genome) 
	//   but we update the direction here respect to the direction of first interval.
	if(!segment.isFwd()){
	    this.graphPos = this.mappedEdge.length() - this.graphPos + 1;
	    this.fwd = !this.fwd;
	}
	}*/

    public void setPairingRead(Read pair){
	this.pairingRead = pair;
    }

    //---------------------------GETTERS----------------------------//

    public String getReadName(){
	return this.readName;
    }

    public Edge getMappedEdge(){
	return this.mappedEdge;
    }
    
    public int getGraphPos(){
	return this.graphPos;
    }
    
    //we always want to return 
    //public int getGraphPos(boolean intervalDirection){
	//if(intervalDirection)
    //return this.graphPos;
	    //else
	    //return this.graphPos + this.length - 1;
    //}
    
    private boolean isFwd(){
	return this.fwd;
    }
    
    public boolean getFwd(){
	return this.fwd;
    }
    
    public Read getPairingRead(){
	return this.pairingRead;
    }
        
    public int getLinPos(){
	return this.linPos;
    }

    public int getLength(){
	return this.length;
    }


    public boolean isFirstRead(){
	return this.firstRead;
    }

    public void setFirstRead(boolean fr){
	this.firstRead = fr;
    }

    public boolean delete(){
	return this.needToDelete;
    }
    
    public void demarkDelete(){
	this.needToDelete = false;
    }

    public void markDelete(){
	this.needToDelete = true;
    }

    private Edge mappedEdge;
    private int graphPos; //distance from 5' side of the forward traversal interval. Must flip for reversal.
    private boolean fwd; /* NEED TO UPDATE THIS DIRECTION so that it's with the respect to the first interval of the mappedEdge */
    private Read pairingRead;
    private int length;
    
    private int linPos;
    
    private String readName;
    private String genome;

    private boolean firstRead;//read1 read2 designation

    private boolean needToDelete;

    private boolean uniq;
    
}


class SimpleReadsLoader{
    
    private ArrayList<Integer> inserts;
    private double[] medMAD;

    SimpleReadsLoader(String samfile, Thread t){
	this.inserts = new ArrayList<Integer>();
	this.medMAD = new double[2];
	this.loadReader(samfile, t);
    }

    
    //check if the pair is --><-- and positive insert size.
    private boolean isCO(String flagStr, int isize){
	if(isize > 0){
	    short testingVal = 0x0030;//0b110000
	    short flag = Short.parseShort(flagStr);
	    //--><--
	    if( (flag & testingVal) == 0x0020)
		return true;
	}
	return false;
    }

    private boolean qcFilter(String mapQual){
	byte mappingQual = Byte.parseByte(mapQual);
	if(mappingQual >= Constants.MIN_MAP_QUAL)
	    return true;
	return false;
    }

    private boolean isMapped(String flagStr){
	short unmapped = 0x0004;
	short flag = Short.parseShort(flagStr);
	if( (flag & unmapped) == 0 )
	    return true;
	return false;
    }

    private boolean isSecondary(String flagStr){
	short notPrimary = 0x0100;
	short flag = Short.parseShort(flagStr);
	if( (flag & notPrimary) == 0 )
	    return false;
	return true;
    }

    //returns next READ from file that is at least qc-passed and mapped
    private Read feed(BufferedReader br){
	String curline = "";
	try{
	    while((curline = br.readLine())!=null){
		if(curline.charAt(0) != '@'){//skip headers
		    String[] tokens = curline.split("\\t");
		    
		    if(Constants.DEBUG4)
			System.err.println("SRL: attempting to load :\t" + tokens[0]);
		    
		    //as long as qcFiltered and mapped to the reference we use the read and return the read obj
		    if(this.qcFilter(tokens[4]) && this.isMapped(tokens[1]) && !isSecondary(tokens[1]) ){
			if(Constants.DEBUG4)
			    System.err.println("\tSRL[QCED & MAPPED]: loading read record" );
			Read r = new Read(tokens);//return as soon as we have next record.
			if(r.getLength() >= Constants.MIN_MAPPINGLEN){
			    //checking if read paid is correctly oriented ---> <--- for iSize distribution inference (added 10/20/15 to infer medMAD internally)
			    int tmpISize = Integer.parseInt(tokens[8]);
			    if(this.isCO(tokens[1], tmpISize))
				this.inserts.add(new Integer(tmpISize));
			    
			    return r;
			}else{
			    if(Constants.DEBUG4)
				System.err.println("\tRL[NO PADD]: skipping this read");
			}
			    
		    }else{//if not mapped or not qc-passed then we discard this read and also pop it's pair if it's in readHash
			if(Constants.DEBUG4)
			    System.err.println("\tRL[NO PADD]: skipping this read" );
		    }
		}
	    }
	    return null; // if we exausted the file return null
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	return null;
    }

    /*
     * Loads each line of reads from a sorted samfile
     * and inserts reads onto the correct location on the linkedList of Edges.
     *
     */
    public void loadReader(String samfile, Thread t){
	Interval curInterval = null;
	Edge curEdge = null;
	BufferedReader br = null;
	String curline = null;
	Read curRead = null;
	
	try{
	    br = new BufferedReader(new FileReader(samfile));
	    curRead = this.feed(br);
	    
	    if(curRead != null){
		
		do{
		    curInterval = t.getIntervalForPositionX(curRead.getLinPos());
		    if(curInterval == null)
			System.err.println(">>>>>READ LIN POS >>>>" + curRead.getLinPos());
		    curRead.updateDepth(curInterval);
		    
		}while( (curRead = this.feed(br)) !=null);
		
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	
	if(!this.writeMedMADFile())
	    System.exit(1);
	    
	t.serializeDepthArray(samfile);
    }

    private double calculateMAD(){
        double [] absDevs = new double[this.inserts.size()];
        for(int i=0; i<this.inserts.size();i++){
            absDevs[i] = Math.abs(this.inserts.get(i).doubleValue() - this.medMAD[0]);
        }
        Arrays.sort(absDevs);
        if(absDevs.length % 2 == 1)
            return absDevs[absDevs.length/2];
        else
            return ( (absDevs[absDevs.length/2] + absDevs[absDevs.length/2 - 1]) / 2.0d );
    }

    private boolean writeMedMADFile(){
	int numData = this.inserts.size();
	if(numData > 0){
	    Collections.sort(this.inserts);
	    if(numData % 2 == 1)
		this.medMAD[0] = this.inserts.get(numData/2)*1.0d;
	    else
		this.medMAD[0] = ( (this.inserts.get(numData/2) + this.inserts.get(numData/2 - 1))*1.0d ) / (2.0d);
	    this.medMAD[1] = this.calculateMAD();
	    
	    BufferedWriter bw = null;
	    try{
		bw = new BufferedWriter(new FileWriter(Constants.PROJNAME + ".medMAD_GRASPER"));
		bw.write(medMAD[0] + "\t" + medMAD[1] + "\n");
		bw.close();
	    }catch(IOException ioe){
		ioe.printStackTrace();
	    }
	}else{
	    this.medMAD[0] = 0.0d;
	    this.medMAD[1] = 0.0d;
	    System.err.println("No corrected oriented paired-end reads to infer medMAD values. Check your reads. System exiting... NO depth file written as a result.");
	    return false;
	}
	return true;
    }    
    
}

//
// load reads from left to rigth from sorted sam file
// in order to pair reads, readHash records the name of the read upon loading. 
// When loading its pair, the name is checked againt the hash to get the proper pairing.
//
class ReadsLoader{

    private HashMap<String, Read> readHash;//simply for pair-name checking to assign pointers to reads in a pair.
    
    public ReadsLoader(){
	this.readHash = new HashMap<String, Read>();
    }

    public ReadsLoader(String samfile, Thread t, boolean updateDepth){
	this();
	this.loadReader(samfile, t, updateDepth);
    }

    private void setPair(Read r1, Read r2){
	r1.setPairingRead(r2);
	r2.setPairingRead(r1);
    }

    //updates the read pair hash
    public Read loadReadRecord(String[] tokens){
	Read curRead = new Read(tokens);
	Read pair = this.readHash.remove(curRead.getReadName());
	if(pair != null){
	    curRead.setFirstRead(false);//this is the second read based on the sorted positions.
	    setPair(pair, curRead);
	}else{
	    curRead.setFirstRead(true);//this is the firest read based on the sorted positions.
	    this.readHash.put(curRead.getReadName(), curRead);
	}
	return curRead;
    }
    
    private boolean qcFilter(String mapQual){
	byte mappingQual = Byte.parseByte(mapQual);
	if(mappingQual >= Constants.MIN_MAP_QUAL)
	    return true;
	return false;
    }

    private boolean isMapped(String flagStr){
	short unmapped = 0x0004;
	short flag = Short.parseShort(flagStr);
	if( (flag & unmapped) == 0 )
	    return true;
	return false;
    }

    private boolean isSecondary(String flagStr){
	short notPrimary = 0x0100;
	short flag = Short.parseShort(flagStr);
	if( (flag & notPrimary) == 0 )
	    return false;
	return true;
    }
    
    //returns next READ from file that is at least qc-passed and mapped
    private Read feedNextRead(BufferedReader br){
	String curline = "";
	try{
	    while((curline = br.readLine())!=null){
		if(curline.charAt(0) != '@'){//skip headers
		    String[] tokens = curline.split("\\t");
		    
		    if(Constants.DEBUG4)
			System.err.println("RL: attempting to load :\t" + tokens[0]);
		    
		    //as long as qcFiltered and mapped to the reference we use the read and return the read obj
		    if(this.qcFilter(tokens[4]) && this.isMapped(tokens[1]) && !isSecondary(tokens[1]) ){
			if(Constants.DEBUG4)
			    System.err.println("\tRL[QCED & MAPPED]: loading read record" );
			return this.loadReadRecord(tokens);//return as soon as we have next record.
		    }else{//if not mapped or not qc-passed then we discard this read and also pop it's pair if it's in readHash
			if(Constants.DEBUG4)
			    System.err.println("\tRL[NO PADD]: skipping this read" );
			if(!isSecondary(tokens[1])){
			    if(Constants.DEBUG4)
				System.err.println("\t\tRL[Primary Hit] : attempt to remove its pair from readHash" );
			    this.readHash.remove(tokens[0]);//pairing implementation for 
			}
		    }

		}
	    }
	    return null; // if we exausted the file return null
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}
	return null;
    }
    
    /*
     * Loads each line of reads from a sorted samfile
     * and inserts reads onto the correct location on the linkedList of Edges.
     *
     */
    public void loadReader(String samfile, Thread t, boolean updateDepth){
	Interval curInterval = null;
	Edge curEdge = null;
	BufferedReader br = null;
	String curline = null;
	Read curRead = null;
	int numReadsAdded = 0;
	String prevReadName = null;
	
	try{
	    br = new BufferedReader(new FileReader(samfile));
	    curRead = this.feedNextRead(br);
	    
	    if(curRead != null){
		
		//for each interval in the thread
		for(int i=0; i<t.getIntervalSequence().size(); i++){
		    curInterval = t.getIntervalSequence().get(i);
		    curEdge = curInterval.getEdge();
		    
		    LinkedList<Read> readList = curEdge.getReads();
		    
		    DirectionalIterator<Read> readIterator = new DirectionalIterator<Read>(readList, curInterval.isFwd());
		    Read fedRead = null;
		    do{ //WITHIN A INTERVAL, we check reads already on the edge to check insertion location
			
			//System.err.println(curRead);

			// 0 : falls within the interval
			// 1 : after the interval
			//-1 : before the interval
			int compareVal = curInterval.compare(curRead.getLinPos());
			
			//curRead belongs to curEdge
			if(compareVal == 0){
			    if(Constants.DEBUG4)
				System.err.println("\tcurRead mapped to the Edge]\t" + curInterval.toString());
			    if(updateDepth)
				curRead.updateDepth(curInterval);
			    else
				curRead.update(curInterval);
			    if(!curRead.isReadRepetitiveOnUniqueEdge()){ //7/7/15 added[heewlee] to prevent mappings with multiple hits with optimal scores on unique edge.
				boolean added = false;
				while(readIterator.hasNext()){
				    fedRead = readIterator.next();
				    if(fedRead == null){
					if(Constants.DEBUG4)
					    System.err.println("Oops: fedRead == null in loadReader(String, Thread)");
				    }
				    //we should insert the curRead as soon as we find a read that is past the position of curRead
				    if(!fedRead.isBeforeOnSameEdgeAtEc(curRead, curInterval.getEC())){
					if(Constants.DEBUG4)
					    System.err.println("Adding curRead right before:\t" + fedRead.toString());
					readIterator.addWhenNotNull(curRead);
					readIterator.previous();
					added = true;
					break;
				    }else
					;//we keep feeding next read from the edge
				}
				if(!added)
				    readIterator.add(curRead);
			    }else{
				curRead.markDelete();
				if(curRead.getPairingRead() != null){
				    curRead.getPairingRead().markDelete();
				    this.readHash.remove(curRead.getReadName());
				}
			    }
			}else if(compareVal == 1) // feed next interval
			    break; // break out of do-while to go to next interval
			else{ // curRead is before the interval. feed next read --> THIS SHOULDN"T HAPPEN!!
			    System.err.println("CRITICAL ERROR IN LOADING READS!!! ORDER MESSED UP :\t" + curRead.getReadName());
			    System.exit(0);
			}
		    }while( (curRead = this.feedNextRead(br)) != null);
		    
		    if(curRead == null)//if not more reads, no need to traverse more
			break; // break out of for to end
		}//end of for
	    
	    }
	    br.close();
	}catch(IOException ioe){
	    ioe.printStackTrace();
	}

	//
	// CODE can be added here to handle what's left in the readHash. Singleton and etc
	//
	this.readHash = null; 

    }


}
