// Interval has reference to the Edge in which the interval belongs to.
//
// Interval also has edge counter value so you know what this interval's EC val is.
//
// Interval has directionality. This direction has nothing to do with the direction 
// of reference strand of ref genome
//
// This directionaltiy only means whether it's forward traversal of the Edge or reverse 
// traversal of the edge. 
//
// The first interval in an intervals array in its Edge, by definition has its
// traversal direction as forward. 
// ==> First traversal of the Edge in the thread is always forward by definition.
//
// EDGE :   A_END ------------------ B_END
// FWD Interval : A_END -3--->>>>>>>> B_END
// REV Interval : A_END <<<<<<<<---- B_END
//
//
// start of interval is A_END if FWD interval
//                      B_END if REV interval

/////////////////////////////////////////////////////////////////////////////////////////////
//    private Edge edge; 
//    private int ec; // 0-based index number: 0->1st occurence, 1->2nd occurence, and so on
//    private boolean fwd;
//    private int start; //5' start of the interval --> 5' is from the reference strand. 
//    private int length; //this length can be different for 2 intervals that are part of single edge.
//    private int threadIndex; //0-based index number for this interval's position on intervalSequence of Thread. --> allows random access index of threading.
/////////////////////////////////////////////////////////////////////////////////////////////
public class Interval implements Comparable<Interval>{

    public int getLinPosFromStartDepthArrIndex(int startIndex){
	return this.getLinearPosition(startIndex*100);
    }

    public int getLinPosFromEndDepthArrIndex(int endIndex){
	return this.getLinearPosition( (endIndex+1) * 100 );
    }
    
    /*public Interval(String line){
	String[] toks = line.split(" ");
	this.fwd = false;
	if(toks[1].charAt(0) == '0')
	    this.fwd = true;
	this.start = Integer.parseInt(toks[2]);
	this.length = Integer.parseInt(toks[3]);
	//toks[4] --> offset : NOT USED. for Al-bruign
	}*/

    public Interval(Edge e, int st, int len, int ec){
	this.edge = e;
	if(st == 0){
	    this.start = st+1;
	    this.length = len;
	}else{
	    this.start = st+2;
	    this.length = len - 1;
	}
	//this.length = length;
	this.ec = ec;
    }
    
    public Interval(Edge e, int dir, int st, int len, int ec){
	this(e,st,len, ec);
	if(dir == 1)
	    this.fwd = false;
	else
	    this.fwd = true;
    }
    
    public Interval(Edge e, boolean direction, int start, int length, int ec){
	this(e,start,length, ec);
	this.fwd = direction;
    }
    
    public void setThreadIndex(int n){
	this.threadIndex = n;
    }

    public int getThreadIndex(){
	return this.threadIndex;
    }
    
    public void setEC(int ec){
	this.ec = ec;
    }

    public Edge getEdge(){
	return this.edge;
    }
    
    public int distanceToNextIntervalFromLinPos(int linpos){
	return this.start + this.length - 1 - linpos;
    }

    public int distanceToPreviousIntervalFromLinPos(int linpos){
	return linpos - this.start;
    }

    public int distanceTo(Interval other){
	if(this.threadIndex == other.getThreadIndex())
	    return 0;
	else if(this.threadIndex < other.threadIndex){
	    return other.getStart() - this.getEnd();
	}else
	    return this.getStart() - other.getEnd();
    }

    public boolean sufficientDistanceToNextInterval(int linPos){
	if(this.length - (linPos - this.start + 1) >= 100)
	    return true;
	return false;
    }

    public boolean sufficientDistanceToPrevInterval(int linPos){
	if( (linPos - this.start + 1) >= 100)
	    return true;
	return false;
    }

    public boolean isSameDirectionIntervalOfSameEdge(Interval other){
	if( this.edge.equals(other.getEdge())
	    && this.fwd == other.fwd )
	    return true;
	else
	    return false;
	    
    }

    public int getEC(){
	return this.ec;
    }
    
    public boolean isFwd(){
	return this.fwd;
    }
    
    public int getEnd(){
	return this.start + this.length() - 1;
    }

    public int getStart(){
	return this.start;
    }
    
    public int length(){
	return this.length;
    }

    /* returns the leftmost position respect to the reference strand */
    public int getLinearPosition(int graphPos){
	return this.start + this.graphPos2IntervalPos(graphPos) - 1;
    }

    /* converts linPos to gp. accounts for directionality as well */
    public int linPos2GraphPos(int linPos){
	int tmp;
	//if length of interval is same as the first interval length (i.e. edge length)
	if(this.length == this.edge.length()){
	    tmp = linPos - this.start + 1;
	}else{//lengths differ!! 
	    //each interval belonging to the same edge can have different(similar) lengths. So here we normalize the length so that it fits within the edge length(length of the first segment)
	    tmp = (int) ( ((linPos - this.start + 1)*1.0d)
			  / (1.0d * this.length)
			  * (this.edge.length() *1.0d)
			  );
	}

	/* for reverse segment, gp has to be flipped. No need to be offset by half-readLength 
	   as linPos used to compute gp is already midpoint initially read direction is loaded
	   according to SAM outputs (that is just respect to the reference genome) 
	   but we update the direction here respect to the direction of first interval.*/
	if(!this.isFwd())
	    return this.edge.length() - tmp + 1;
	return tmp;
    }

    public int graphPos2IntervalPos(int graphPos){
	int tmp = (int) (((graphPos*1.0d)/(this.edge.length()*1.0d))*(this.length*1.0d));
	if(this.fwd)
	    return tmp;
	else
	    return this.length - tmp + 1;
    }
    
    
    public int linPos2RelativePosOnInterval(int linPos){
	return linPos - this.start + 1;
    }

    public Interval getNextInterval(){
	return this.edge.nthInterval(this.ec+1);
    }
    
    public boolean isLastInterval(){
	if( this.ec == (this.edge.getMultiplicity()-1) )
	    return true;
	return false;
    }

    public boolean isSameInterval(Interval other){
	if( this.threadIndex == other.getThreadIndex() )
	    return true;
	return false;
    }

    public boolean coverLinPos(int linPos){
	if(linPos >= this.start 
	   && linPos < (this.start + this.length) )
	    return true;
	return false;
    }
        
    /* returns -1,0,1 if x is smaller/within/bigger*/
    public int compare(int x){
	if(x >= this.start){
	    if(x <= (this.start+this.length-1))
		return 0;
	    return 1;
	}else
	    return -1;
    }
    
    //ascending order
    public int compareTo(Interval other){
	return this.start - other.start;
    }

    public String toString(){
	return "edge:" + this.edge.getId() + "\t" + "ec:" + this.ec + "\t" + "fwd:"
	    + this.fwd + "\t" + "start:" + this.start + "\t" + "length:" + this.length;
    }

    //////////////////////////////////
    private Edge edge; 
    private int ec; // 0-based index number: 0->1st occurence, 1->2nd occurence, and so on
    private boolean fwd;
    private int start; //5' start of the interval
    private int length; //this length can be different for 2 intervals that are part of single edge.
    private int threadIndex; //0-based index number for this interval's position on intervalSequence of Thread. --> allows random access index of threading.
    //////////////////////////////////
}

