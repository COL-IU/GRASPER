import java.util.*;

// Edge has unique STRING identifier called id
// Edge has array of intervals where the length of the array is equal to the multiplicty of the edge.
// Edge can have both directional traversal as long as intervals in intervals-array has both direction
// Edge can be thought of as having identical direction as the reference strand of ref genome. 
// 
// Edge has a set of reads and a set of clusters --> data structure is set as arrayList as of now but linked list
// can be considered.
//
///////CLASS FIELDS//////////////////
/////////////////////////////////////
//    private String id;
//    private Interval[] intervals;
//    private LinkedList<Read> reads;
//    private LinkedList<Cluster> clusters;
//    private int[] depthArrFwd;
//    private int[] depthArrRev;
////////////////////////////////////
//
public class Edge{

    public Edge(String i, /*int l,*/ int m){
	this.id = i;
	//this.multiplicity = m;
	this.intervals = new Interval[m];
	this.reads = new LinkedList<Read>();//new ArrayList<Read>();
	this.clusters = new LinkedList<HalfCluster>();
	this.depthArrFwd = null;
	this.depthArrRev = null;
    }

    public int[] getDepthArrFwd(){
	return this.depthArrFwd;
    }
    
    public void setDepthArrFwd(int[] f){
	this.depthArrFwd = f;
    }
    
    public int[] getDepthArrRev(){
	return this.depthArrRev;
    }
    
    public void setDepthArrRev(int[] r){
	this.depthArrRev = r;
    }
    

    //@Override
    public int hashCode(){
	return this.id.hashCode();
    }
    //@Override
    public boolean equals(Object obj){
	if(this == obj)
	    return true;
	if(obj == null)
	    return false;
	if(getClass() != obj.getClass())
	    return false;
	Edge other = (Edge) obj;
	return this.id.equals(other.getId());
    }
    
    public String toString(){
	return "Edge:" + id + "Length:" + this.length() +"|numIntervals:" + intervals.length+ "|depthArrLength:" + depthArrFwd.length + ":" + depthArrRev.length;
    }

    /* i : id of the edge 
     * m : multiplicity
     * tokens: tab-tokens of the intervalLine in the thread file
     */
    public Edge(String i, /*int l,*/int m, String[] tokens, int genomeLength){
	this(i, m);
	
	if(!loadIntervals(tokens, genomeLength, m)){
	    System.err.print("INCONSISTENT DATA [LOADING INTERVALS FOR AN EDGE] :\n -->");
	    for(int k=0;k<tokens.length;k++)
		System.err.print("\t" + tokens[k]);
	    System.err.println("===========>>> EXITING PROGRAM ");
	    System.exit(-1);
	}
	
	this.initDepthArrs();
    }

    //public boolean isBefore(Edge other){
    //	;
    //}

    public boolean isUniqueEdge(){
	if(this.intervals.length > 1)
	    return false;
	return true;
    }
    
    private boolean initDepthArrs(){
	if(this.intervals[0] != null){
	    if(this.length()%Constants.DEPTH_WIN_SIZE == 0){
		this.depthArrFwd = new int[this.length()/Constants.DEPTH_WIN_SIZE];
		this.depthArrRev = new int[this.length()/Constants.DEPTH_WIN_SIZE];
	    }else{
		this.depthArrFwd = new int[this.length()/Constants.DEPTH_WIN_SIZE + 1];
		this.depthArrRev = new int[this.length()/Constants.DEPTH_WIN_SIZE + 1];
	    }
	    return true;
	}
	return false;
    }

    public void incrementDepthAt(int gp, boolean directionRespectToEdge){
	int index = this.gpToDepthArrIndex(gp);//(gp-1) / Constants.DEPTH_WIN_SIZE;
	try{
	    if(directionRespectToEdge)
		this.depthArrFwd[index]++;
	    else
		this.depthArrRev[index]++;
	}catch(Exception e){
	    System.err.println(this.toString());
	    System.err.println("GP:" + gp + "|index:" + (gp/Constants.DEPTH_WIN_SIZE));
	    e.printStackTrace();
	    
	}
    }

    public int getDepthAtIndex(int depthIndex){
	return this.depthArrFwd[depthIndex] + this.depthArrRev[depthIndex];
    }
    
    private int gpToDepthArrIndex(int gp){
	int index = (gp - 1) / Constants.DEPTH_WIN_SIZE;
	if(index >= this.depthArrFwd.length)
	    return -1;
	else
	    return index;
    }
    
    public int getTotalDepth(int startgp, int endgp){
	int sIndex = -1;
	int eIndex = -1;
	if(startgp > endgp){
	    sIndex = this.gpToDepthArrIndex(endgp);
	    eIndex = this.gpToDepthArrIndex(startgp);
	}else{
	    sIndex = this.gpToDepthArrIndex(startgp);
	    eIndex = this.gpToDepthArrIndex(endgp);
	}
	
	int total = 0;
	for(int i=sIndex; i<=eIndex; i++){
	    total += this.getDepthAtIndex(i);
	}
	return total;
    }
    

    /*
     * returns:
     *
     * 1  : if edge has insufficient coverages, regardless of whether Edge is repetative or not
     * 0  : if edge has coverages and edge has to be non-repetative
     * -1 : if edge has coverages and edge has to be repetative
     * -2 : gp1 and gp2 too close (no depthArray between two points)
     */
    public int isPartialEdgeLackingCoverage(int gp1, int gp2){
	int sIndex = 0;
	int eIndex = 0;
	if(gp2 < gp1){
	    sIndex = this.gpToDepthArrIndex(gp2);
	    eIndex = this.gpToDepthArrIndex(gp1);
	}else{
	    sIndex = this.gpToDepthArrIndex(gp1);
	    eIndex = this.gpToDepthArrIndex(gp2);
	}

	//we are checking if there is at least one index between sIndex and eIndex.
	//we require at least one cell to validate the depth.
	if( (eIndex - sIndex) > 1 ){
	    	    
	    int insuffDepthBinsCount = 0;
	    
	    for(int i=sIndex+1; i<eIndex; i++){
		if(this.getDepthAtIndex(i) < Constants.MIN_DEPTH)
		    insuffDepthBinsCount++;
	    }
	    double lackingBinsRatio = (insuffDepthBinsCount*1.0d) / ((eIndex - sIndex - 1)*1.0d);
	    if(lackingBinsRatio > Constants.NO_COVERAGE_RATIO)
		return 1;
	    else if(this.isRepetitive())
		return -1;
	    else
		return 0;
	}else{
	    return -2;
	}
    }

    public int computeSegmentLength(int gp, boolean intervalDirection, boolean uptoGP){
	if(intervalDirection != uptoGP)
	    return (this.length() - gp + 1);
	else
	    return gp;
    }
    
    /*
     * returns an array size of 3. 0: lackingVal, 1: insuffDepthBinsCount, 2: totalBinsCount
     *
     * 1  : if edge has insufficient coverages, regardless of whether Edge is repetative or not
     * 0  : if edge has coverages and edge has to be non-repetative
     * -1 : if edge has coverages and edge has to be repetative
     * 
     */
    public int[] isPartialEdgeLackingCoverage(int gp, boolean intervalDirection, boolean uptoGP){
	int[] result = new int[3];

	int insuffDepthBinsCount = 0;
	int sIndex = this.gpToDepthArrIndex(gp);
	
	//  IF interval dir is same as edge direction and from gp to end
	//OR  IF interval dir REV, and from beginning of interval to GP.  --> from gp to end on EDGE
	if(intervalDirection != uptoGP){ // depth from left side of the edge :|||||||||(gp)--------
	    for(int i=sIndex+1; i<this.depthArrFwd.length;i++){
		if(this.getDepthAtIndex(i) < Constants.MIN_DEPTH)
		    insuffDepthBinsCount++;
	    }
	    result[1] = insuffDepthBinsCount;
	    result[2] = this.depthArrFwd.length - sIndex -1;
	    
	    double lackingBinsRatio = (result[1]*1.0d) / ((result[2])*1.0d);
	    if(lackingBinsRatio > Constants.NO_COVERAGE_RATIO){
		result[0] = 1;
		return result;
	    }
	}else{ // // depth from right side of the edge :--------(gp)||||||||
	    for(int i=sIndex-1; i>=0;i--){
		if(this.getDepthAtIndex(i) < Constants.MIN_DEPTH)
		    insuffDepthBinsCount++;
	    }
	    result[1] = insuffDepthBinsCount;
	    result[2] = sIndex;
	    double lackingBinsRatio = (insuffDepthBinsCount*1.0d) / (sIndex*1.0d);
	    if(lackingBinsRatio > Constants.NO_COVERAGE_RATIO){
		result[0] = 1;
		return result;
	    }
	}

	
	if(this.isRepetitive())
	    result[0] = -1;
	else
	    result[0] = 0;
	return result;
    }
   
    public boolean isRepetitive(){
	if(this.intervals.length > 1)
	    return true;
	return false;
    }

    /*
     * returns:
     *
     * 1  : if edge has insufficient coverages, regardless of whether Edge is repetative or not
     * 0  : if edge has coverages and edge has to be non-repetative
     * -1 : if edge has coverages and edge has to be repetative
     * 
     */
    public int[] isEdgeLackingCoverage(){
	int[] result = new int[3];

	int insuffDepthBinsCount = 0;
	for(int i=0;i<this.depthArrFwd.length;i++){
	    if(this.getDepthAtIndex(i)  < Constants.MIN_DEPTH)
		insuffDepthBinsCount++;
	}
	result[1] = insuffDepthBinsCount;
	result[2] = this.depthArrFwd.length;
	
	double lackingBinsRatio = (insuffDepthBinsCount*1.0d) / (this.depthArrFwd.length*1.0d);
	if(lackingBinsRatio > Constants.NO_COVERAGE_RATIO)
	    result[0] = 1;
	else if(this.isRepetitive())
	    result[0] = -1;
	else
	    result[0] = 0;
	return result;
    }
    
    public double ratioOfInsufficientDepthBins(int startgp, int endgp){
	if(this.getMultiplicity() == 1)
	    return 0.0d;
	
	int sIndex = -1;
	int eIndex = -1;
	if(startgp > endgp){
	    sIndex = this.gpToDepthArrIndex(endgp);
	    eIndex = this.gpToDepthArrIndex(startgp);
	}else{
	    sIndex = this.gpToDepthArrIndex(startgp);
	    eIndex = this.gpToDepthArrIndex(endgp);
	}
	
	//int total = 0;
	int insuffDepthBinsCount = 0;
	for(int i=sIndex; i<=eIndex; i++){
	    if(this.getDepthAtIndex(i) < Constants.MIN_DEPTH)
		insuffDepthBinsCount++;
	}
	return (insuffDepthBinsCount*1.0d)/((eIndex-sIndex+1)*1.0d);
    }

    
    private boolean loadIntervals(String[] tokens, int genomeLength, int m){
	/* multiplicity must be equal to interval contained in the threadFile */
	if( (tokens.length - 4)/3 != m ) 
	    return false; //inconsistent data
	
	
	/* 
	 *  field 4i   : direction --> direction 0:forward, 1:reverse
	 *  field 4i+1 : direction wise start --> if reverse interval, genomeLength - start - length will get fwd position
	 *  field 4i+2 : length
	 */
	for(int i=4; i<tokens.length; i=i+3){
	    this.intervals[(i-4)/3] = this.loadInterval( Integer.parseInt(tokens[i]) ,
							 Integer.parseInt(tokens[i+1]) ,
							 Integer.parseInt(tokens[i+2]) ,
							 genomeLength ,
							 (i-4)/3 );
	}
	
	/* need to sort the intervals in the order of fowrard thread: 
	 * that is the direction of reference strand 5' -> 3'*/
	this.sortIntervals(); 
	
	return true;
    }
    

    
    private Interval loadInterval(int direction, int start, int length, int genomeLength, int ec){
	/* if direction is reverse */
	if(direction == 1){
	    start = genomeLength - start - length; /* getting the coordinate for the other end */
	    return new Interval(this, false, start, length, ec);
	}else
	    return new Interval(this, true, start, length, ec); 
    }
    
    /* sorts the intervals arrays */
    public void sortIntervals(){
	Arrays.sort(this.intervals);
	this.updateECs();//after sorting EC (Edge Counter) needs to be updated.s
    }
    
    /* EC val for Interval in Intervals array is identical to the array index */
    private void updateECs(){
	for(int i=0;i<this.intervals.length; i++)
	    this.intervals[i].setEC(i);
    }
    
    // graphPos is the graphPos belonging to THIS edge. We want to compute the shortest distance to read R.
    /*public int getDistanceFromGraphPosTo(int graphPos, Read r){
	int thisMultiplicity = this.getMultiplicity();
	int otherMultiplicity = r.getMappedEdge().getMultiplicity();
	
	if(thisMultiplicity <= otherMultiplicity){
	    
	}else{
	
	}
	}*/
    
    /// NNED TO CHEKC THIS OVER
    public int getDistanceFromGraphPosTo(int graphPos, Read r, int maxDistance){
	return Edge.getDistanceBetween2GPs(this, graphPos, r.getMappedEdge(), r.getGraphPos(), maxDistance);
    }

    
    public int[] getLinearPositions(int graphPos){
	int[] linPoss = new int[this.intervals.length];
	
	for(int i=0;i<this.intervals.length;i++)
	    linPoss[i] = intervals[i].getLinearPosition(graphPos);
	
	return linPoss;
    }

    public String possibleLinearPositionPairs(Read minGPRead, Read maxGPRead, boolean direction){
	//public String possibleLinearPositionPairs(int minGP, int maxGP, boolean direction){
	int minGP = minGPRead.getGraphPos();
	int maxGP = maxGPRead.getGraphPos();
	StringBuffer buffer = new StringBuffer();
	for(int i=0;i<this.intervals.length;i++){
	    int tmp1 = intervals[i].getLinearPosition(minGP);
	    int tmp2 = intervals[i].getLinearPosition(maxGP);
	    
	    if(tmp1<tmp2)
		buffer.append((tmp1 - minGPRead.getLength()/2) + "\t" + (tmp2 + maxGPRead.getLength()/2) + "\t");
	    else
		buffer.append( (tmp2 - maxGPRead.getLength()/2) + "\t" + (tmp1 + minGPRead.getLength()/2) + "\t");
	    
	    buffer.append(( (direction == intervals[i].isFwd()) ? "--->" : "<---") + "\n");
		
	}
	return buffer.toString();
    }
    
    public DirectionalRange getFirstLinearPositionPairsAsDirectionalRange(Read minGPRead, Read maxGPRead, boolean direction){
	return getNthDirectionalRange(minGPRead, maxGPRead, 0);
    }
    
    /* direction is already respect to the genome at n-th interval */
    //Taking the direction from min/maxGPRead
    public DirectionalRange getNthDirectionalRange(Read minGPRead, Read maxGPRead, int n){
	boolean direction1 = minGPRead.isFwdAtEc(n);
	boolean direction2 = maxGPRead.isFwdAtEc(n);
	boolean direction = direction1;
	if(direction1 != direction2){
	    System.err.println("##############################################################\n\t\t\t\tOOPS DIRECTIONS DIFFER in edgeCLUST???\n############################################################");
	    System.exit(1);
	}
	
	int minGP = minGPRead.getGraphPos();
	int maxGP = maxGPRead.getGraphPos();
	int tmp1 = intervals[n].getLinearPosition(minGP);
	int tmp2 = intervals[n].getLinearPosition(maxGP);
	if(tmp1<tmp2){
	    return new DirectionalRange( (tmp1 - minGPRead.getLength()/2) 
					 , (tmp2 + maxGPRead.getLength()/2)
					 , direction
					 , intervals[n].getThreadIndex());//(direction == intervals[n].isFwd()) );
	}else
	    return new DirectionalRange( (tmp2 - maxGPRead.getLength()/2)
					 , (tmp1 + minGPRead.getLength()/2)
					 , direction
					 , intervals[n].getThreadIndex());//(direction == intervals[n].isFwd()) );
    }


    public String getFirstLinearPositionPairs(Read minGPRead, Read maxGPRead, boolean direction){
	//public String getFirstLinearPositionPairs(int minGP, int maxGP, boolean direction){
	int minGP = minGPRead.getGraphPos();
	int maxGP = maxGPRead.getGraphPos();
	StringBuffer buffer = new StringBuffer("\t" + this.getMultiplicity());
	int tmp1 = intervals[0].getLinearPosition(minGP);
	int tmp2 = intervals[0].getLinearPosition(maxGP);

	if(tmp1<tmp2)
	    buffer.append("\t" + (tmp1 - minGPRead.getLength()/2)  + "\t" + (tmp2 + maxGPRead.getLength()/2));
	else
	    buffer.append("\t" + (tmp2 - maxGPRead.getLength()/2)  + "\t" + (tmp1 + minGPRead.getLength()/2));
	buffer.append(( (direction == intervals[0].isFwd()) ? "\t--->" : "\t<---") );
	buffer.append("\tEDGE:" + this.getId());
	return buffer.toString();
    }
    // computes mindistance between two points exaustively
    // returns a negative value (-1) if dist was not less than Constants.DM
    public static int getDistanceBetween2GPs(Edge e1, int gp1, Edge e2, int gp2, int maxDistance){
	int j=0;
	int minDist = maxDistance + 1; //since we want to get distance upto DM val.
	boolean matchFound = false;
	for(int i=0; i< e1.getIntervals().length && j<e2.getIntervals().length; i++){
	    int e1_curmin = e1.getIntervals()[i].getStart();
	    int e2_curmin = e2.getIntervals()[j].getStart();
	    if(e2_curmin > e1_curmin){
		if(i < (e1.getIntervals().length - 1) ){
		    for(i=i+1; i<e1.getIntervals().length;i++){
			int e1_innerCurmin = e1.getIntervals()[i].getStart();
			if(e2_curmin < e1_innerCurmin){
			    //i = i -1;
			    break;
			}
		    }
		    i= i - 1;
		    e1_curmin = e1.getIntervals()[i].getStart();
		}
		int dist = e2.getIntervals()[j].getLinearPosition(gp2) 
		    - e1.getIntervals()[i].getLinearPosition(gp1);//no longer needing readLen offset as we do midpoint comparison. // + Constants.READLEN;
		if(dist < minDist){
		    minDist = dist;
		    matchFound = true;
		}
	    }else{ // if e2_curmin comes before --> we need to crack up e2s
		
		if(j < (e2.getIntervals().length - 1) ){
		    for(j=j+1; j<e2.getIntervals().length; j++){
			int e2_innerCurmin = e2.getIntervals()[j].getStart();
			if(e1_curmin < e2_innerCurmin){
			    //j = j - 1;
			    break;
			}
		    }
		    j=j-1;
		    e2_curmin = e2.getIntervals()[j].getStart();
		}
		int dist = e1.getIntervals()[j].getLinearPosition(gp1) 
		    - e2.getIntervals()[i].getLinearPosition(gp2);//no longer needing readLen offset as we do midpoint comparison. // + Constants.READLEN;
		if(dist < minDist){
		    minDist = dist;
		    matchFound = true;
		}else{
		    i--;
		    j++;
		}
	    }
	}
	if(matchFound)
	    return minDist;
	else
	    return -1;
    }

    //boolean version to shorten the time.
    public static boolean isDistanceBetween2GPsWithinDM(Edge e1, int gp1, Edge e2, int gp2, int maxDistance){
	if(e1.equals(e2))
	    return ( (Math.abs(gp2-gp1+1) <= maxDistance) ? true : false );
	
	int j=0;
	for(int i=0; i< e1.getIntervals().length && j<e2.getIntervals().length; i++){
	    int e1_curmin = e1.getIntervals()[i].getStart();
	    int e2_curmin = e2.getIntervals()[j].getStart();
	    if(e2_curmin > e1_curmin){
		if(i < (e1.getIntervals().length - 1) ){
		    for(i=i+1; i<e1.getIntervals().length;i++){
			int e1_innerCurmin = e1.getIntervals()[i].getStart();
			if(e2_curmin < e1_innerCurmin){
			    //i = i - 1;
			    break;
			}
		    }
		    i=i-1;
		    e1_curmin = e1.getIntervals()[i].getStart();
		}
		int dist = e2.getIntervals()[j].getLinearPosition(gp2) 
		    - e1.getIntervals()[i].getLinearPosition(gp1);//no longer needing readLen offset as we do midpoint comparison. // + Constants.READLEN;
		if(dist <= maxDistance)
		    return true;
	    }else{ // if e2_curmin comes before --> we need to crack up e2s
		if(j < (e2.getIntervals().length - 1) ){
		    for(j=j+1; j<e2.getIntervals().length; j++){
			int e2_innerCurmin = e2.getIntervals()[j].getStart();
			if(e1_curmin < e2_innerCurmin){
			    //j = j - 1;
			    break;
			}
		    }
		    j=j-1;
		    e2_curmin = e2.getIntervals()[j].getStart();
		}
		int dist = e1.getIntervals()[i].getLinearPosition(gp1) 
		    - e2.getIntervals()[j].getLinearPosition(gp2);//no longer needing readLen offset as we do midpoint comparison. // + Constants.READLEN;
		if(dist <= maxDistance)
		    return true;
		else{
		    i--;
		    j++;
		}
	    }
	}
	return false;
    }

    //
    // QUERYING POSTION: toPosition (linear position)
    // TARGET POSITION: posOnInterval (graphPos)
    //
    // ------intervalA---posOnInterval----->   ----------toPosition----->   ------intervalA+1----posOnInterval----->
    //
    // whichever betwween interval A and (A+1), that is closer to 'toPosition' is returned.
    public Interval getClosestInterval(int posOnInterval, int toPosition){
	int max = this.intervals.length;
	int min = 0;
	int upper = 0;
	int prevMax = 0;
	while(max >= min){
	    int curIndex = (min+max)/2;
	    int compVal = this.intervals[curIndex].compare(toPosition);
	    //prevMax = max;
	    if(compVal < 0)//position is before the interval.
		max = curIndex - 1;
	    else if (compVal > 0)
		min = curIndex + 1;
	    else
		return this.intervals[curIndex];
	}
	//two intervals at indices max & max+1
	//return the closer one of the two.
	return this.getCloserInterval(posOnInterval, max, max+1, toPosition);
    }

    
    public int getDistanceToClosestInterval(int posOnInterval, int toPosition){
	int max = this.intervals.length;
	int min = 0;
	int upper = 0;
	int prevMax = 0;
	while(max >= min){
	    int curIndex = (min+max)/2;
	    int compVal = this.intervals[curIndex].compare(toPosition);
	    //prevMax = max;
	    if(compVal < 0)//position is before the interval.
		max = curIndex - 1;
	    else if (compVal > 0)
		min = curIndex + 1;
	    else
		return Math.abs(toPosition - this.intervals[curIndex].getLinearPosition(posOnInterval)) + 1;
	}
	//two intervals at indices max & max+1
	//return the closer one of the two.
	return this.getDistanceToCloserInterval(posOnInterval, max, max+1, toPosition);
    }

    private int getDistanceToCloserInterval(int graphPos, int ecLeft, int ecRight, int toPosition){
	int dFromLeftInterval = Math.abs(toPosition - this.intervals[ecLeft].getLinearPosition(graphPos)) + 1;
	int dFromRightInterval = Math.abs(toPosition - this.intervals[ecRight].getLinearPosition(graphPos)) + 1;
	if(dFromLeftInterval <= dFromRightInterval)
	    return dFromLeftInterval;
	else
	    return dFromRightInterval;
    }

    private Interval getCloserInterval(int graphPos, int ecLeft, int ecRight, int toPosition){
	int dFromLeftInterval = Math.abs(toPosition - this.intervals[ecLeft].getLinearPosition(graphPos)) + 1;
	int dFromRightInterval = Math.abs(toPosition - this.intervals[ecRight].getLinearPosition(graphPos)) + 1;
	if(dFromLeftInterval <= dFromRightInterval)
	    return this.intervals[ecLeft];
	else
	    return this.intervals[ecRight];
    }

    //DONT NEED TO THIS HANDELD BY coordinate calculation methods done in Interval class.
    /*
    private int getDistanceToCloserInterval(int posOnInterval, int ecLeft, int ecRight, int toPosition){
	int dTo3PrimeLoc = this.getDistanceTo3PrimePosition(posOnInterval, ecLeft, toPosition);
	int dTo5PrimeLoc = this.getDistanceTo5PrimePosition(posOnInterval, ecRight, toPosition);
	if( dTo3PrimeLoc < dTo5PrimeLoc )
	    return dTo3PrimeLoc;
	else
	    return dTo5PrimeLoc;
    }

    private Interval getCloserInterval(int posOnInterval, int ecLeft, int ecRight, int toPosition){
	if( this.getDistanceTo3PrimePosition(posOnInterval, ecLeft, toPosition) 
	    < this.getDistanceTo5PrimePosition(posOnInterval, ecRight, toPosition) )
	    return this.intervals[ecLeft];
	else
	    return this.intervals[ecRight];
    }

    
    //fwd : --------------<posOnInterval>----->    toPosition    :   toPosition - (this.start + posOnInterval)
    //rev : <-------------<posOnInterval>------    toPosition    :   toPosition - (this.start + length - posOnInterval)   
    private int getDistanceTo3PrimePosition(int posOnInterval, int ec, int toPosition){
	if(this.intervals[ec].isFwd())
	    return toPosition - (this.intervals[ec].getStart() + posOnInterval);
	else
	    return toPosition - (this.intervals[ec].getStart() + this.intervals[ec].length() - posOnInterval);
    }

    //fwd : toPosition    [start]------------------<posOnInterval>----[start + length]->
    //rev : toPosition    <-[start + length]]----------------<posOnInterval>-----[start]
    private int getDistanceTo5PrimePosition(int posOnInterval, int ec, int toPosition){
	if(this.intervals[ec].isFwd())
	    return (this.intervals[ec].getStart() + posOnInterval) - toPosition;
	else
	    return (this.intervals[ec].getStart() + this.intervals[ec].length() - posOnInterval) - toPosition;
    }
    */
    
    /* length of Edge is always equals to the first occurence of the Edge in the thread --> intervals[0] */
    public int length(){
	return this.intervals[0].length();
    }

    public String getId(){
	return this.id;
    }

    public int getMultiplicity(){
	return this.intervals.length;
    }
    
    public Interval[] getIntervals(){
	return this.intervals;
    }
    
    public Interval nthInterval(int n){
	if(n<0 && n>= this.intervals.length)
	    return null;
	return this.intervals[n];
    }
    
    public void addRead(Read r){
	this.reads.add(r);
    }

    
    public LinkedList<Read> getReads(){
	return this.reads;
    }

    public LinkedList<HalfCluster> getClusters(){
	return this.clusters;
    }

    public boolean areReadsSameDirectionClosestIntervals(Edge other, Read tr, Read or){
	boolean withinDM = false;
	int i=0;
	int j=0;
	while(i<this.getMultiplicity() && j<other.getMultiplicity()){
	    //thisPartInterval is before otherPartInterval
	    if(this.nthInterval(i).getEnd() < other.nthInterval(j).getEnd()){
		int diff = other.nthInterval(j).getStart() - this.nthInterval(i).getEnd();
		if(diff < Constants.MAX_D_MIDPOINTS){
		    withinDM = true;
		    if(tr.isFwdAtEc(i) == or.isFwdAtEc(j))
			return true;//sameDirection = true;
		}
		i++;
	    }
	    //otherPartInterval is before thisPartInterval
	    else{ 
		int diff = this.nthInterval(i).getStart() - other.nthInterval(j).getEnd();
		if(diff < Constants.MAX_D_MIDPOINTS){
		    withinDM = true;
		    if(tr.isFwdAtEc(i) == or.isFwdAtEc(j))
			return true;
		}
		j++;
	    }
	}
	
	if(withinDM){
	    return false;
	}
	return true;
    }



    ///////CLASS FIELDS//////////////////
    /////////////////////////////////////
    private String id;
    private Interval[] intervals;
    //private ArrayList<Read> reads;
    //private ArrayList<Cluster> clusters;
    
    private LinkedList<Read> reads;
    private LinkedList<HalfCluster> clusters;
    //private int length;
    //private int multiplicity;
    //private LinkedList<String> reads;
    
    private int[] depthArrFwd; //length of this array is depending on the edge length (length of the first interval belonging to the edge
    private int[] depthArrRev;
    ////////////////////////////////////
}

