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

public class DirectionalRange{

    private int min;
    private int max;
    private boolean fwd;
    private int tipIntervalIndex;  /* this is the array index to access the tip Interval in HalfCluster: index used here is the threadIndex*/
    private int stemIntervalIndex; /* this is the array index to access the stem Interval in HalfCluster */
    //private HalfCluster cl;
    
    //private DirectionalRange pair;
    
    public DirectionalRange(int m, int M){
	this.min = m;
	this.max = M;
    }

    public DirectionalRange(int m, int M, boolean f, int i){
	this(m, M, f, i, i);
    }

    public DirectionalRange(int m, int M, boolean f, int ti, int si){
	this.min = m;
	this.max = M;
	this.fwd = f;
	this.tipIntervalIndex = ti;
	this.stemIntervalIndex = si;
    }


    public int getInversionType(){
	if(this.min == this.max)
	    return this.min;
	else{
	    System.exit(0);
	    return -1;
	}
    }
    
    //boundaries int[6] 
    //        deletion boundary  (int[0], int[1])
    //        5' repeat boundary (int[2], int[3])
    //        3' repeat boundary (int[4], int[5])
    public String toDeletionRangeRepeat(int[] boundaries){
	return "DELETION\tDONPFR\t-\t-\t-\t-\t" + this.min + "\t" + this.max + "\t5'-REPEAT(" +boundaries[2] + "," + boundaries[3] + ")\t3'-REPEAT(" + boundaries[4] + "," + boundaries[5] + ")";
    }
    
    public String toDeletionRange(){
	return "DELETION\tDONP\t-\t-\t-\t-\t" + this.min + "\t" + this.max;
    }
    
    //
    // returns -1 : if oDR is before
    // returns  0 : if oDR is contained 
    // returns  1 : if oDR is after
    public int contains(DirectionalRange oDR){
	if(oDR.getMin() < (this.min - Constants.DEPTH_WIN_SIZE) )
	    return -1;
	else if( (this.min - Constants.DEPTH_WIN_SIZE) <= oDR.min && (this.max + Constants.DEPTH_WIN_SIZE) >= oDR.max)
	    return 0;
	else
	    return 1;
	//return false;
    }

    public String toString(){
	return this.min + "\t" + this.max + "\t" + (this.fwd ? "---->" : "<----");
    }
    
    public String toDeletionStringNoCluster(){
	return this.min + "\t" + this.max;
    }


    public int getMin(){
	return this.min;
    }
    
    public int getMax(){
	return this.max;
    }
    
    public boolean isFwd(){
	return this.fwd;
    }

    public int getTipIntervalIndex(){
	return this.tipIntervalIndex;
    }
    
    public int getStemIntervalIndex(){
	return this.stemIntervalIndex;
    }

    public void setTipIntervalIndex(int n){
	this.tipIntervalIndex = n;
    }
    
    public void setStemIntervalIndex(int n){
	this.stemIntervalIndex = n;
    }

    // --->    --->  or  <---     <---
    //return 0 , 1 , 2 , 3
    public int isFollower(DirectionalRange pair){
	
	if(this.fwd && pair.isFwd()){
	    //---this--->    ---pair--->
	    if(this.max < pair.getMin())
		return 0;
	    //---pair--->    ---this---->
	    else if(pair.getMax() < this.min)
		return 3;
	}else if(!this.fwd && !pair.isFwd()){
	    //<---this---     <---pair---
	    if(this.max < pair.getMin())
		return 1;
	    //<---pair---     <---this---
	    else if(pair.getMax() < this.min)
		return 2;
	}
	return -1;
    }

    /*
     *  0: [NO, TP, ++] --this--> --pair-->
     *  1: [NO, TP, --] <--this-- <--pair--
     *  2: [NO, PT, --] <--pair-- <--this--
     *  3: [NO, PT, ++] --pair--> --this-->
     *
     *  4: [NO, TP, +-] --this--> <--pair--
     *  5: [NO, TP, -+] <--this-- --pair-->
     *  6: [NO, PT, +-] --pair--> <--this--
     *  7: [NO, PT, -+] <--pair-- --this-->
     
     *  8: [OV, TP, ++] -=this,pair=>>
     *  9: [OV, TP, --] <<=this,pair=- 
     * 10: [OV, PT, --] <<=pair,this=-
     * 11: [OV, PT, ++] -=pair,this=>>
     *
     * 12: [OV, TP, +-] -<P=this,pair=T>-
     * 13: [OV, TP, -+] -<T=this,pair=P>-
     * 14: [OV, PT, +-] -<T=pair,this=P>-
     * 15: [OV, PT, -+] -<P=pair,this=T>-
     *
     * -1: OTHERWISE
     */
    public int getRelativeRelationType(DirectionalRange pair){
	// -> ->
	if(this.fwd && pair.isFwd()){
	    /* NON-OVERLAPPING */

	    //---this--->   ---pair--->
	    if(this.max < pair.getMin()){
		return 0;
	    }
	    //---pair--->   ---this--->
	    else if(pair.getMax() < this.min)
		return 3;

	    /* OVERLAPPING */
	    //  this-pair  >>
	    else if(pair.getMax() > this.max)
		return 8;
	    // pair-this >>
	    else if(this.max >= pair.getMax())
		return 11;
	}
	// <- <-
	else if(!this.fwd && !pair.isFwd()){
	    /* NON-OVERLAPPING */

	    //<---this---     <---pair---
	    if(this.max < pair.getMin())
		return 1;
	    //<---pair---     <---this---
	    else if(pair.getMax() < this.min)
		return 2;

	    /* OVERLAPPING */
	    //  this-pair  <<
	    else if(pair.getMin() > this.min)
		return 9;
	    // pair-this <<
	    else if(this.min >= pair.getMin())
		return 10;
	    
	}
	// -> <- or <- ->
	else if(this.fwd && !pair.isFwd()){
	    /* NON-OVERLAPPING */

	    // ---this--->   <----pair---
	    if(this.max < pair.getMin())
		return 4;
	    // <--pair---    ----this--->
	    else if(pair.getMax() < this.min)
		return 7;

	    /* OVERLAPPING */
	    // --this--->
	    //  <---pair---
	    if(this.min < pair.getMin())
		return 12;
	    //  --this--->
	    //<--pair---
	    else if(this.min >= pair.getMin())
		return 15;
	}
	// -> <- or <- ->
	else if(!this.fwd && pair.isFwd()){
	    /* NON-OVERLAPPING */
	    // <---this---   ----pair--->
	    if(this.max < pair.getMin())
		return 5;
	    // --pair--->    <----this---
	    else if(pair.getMax() < this.min)
		return 6;

	    /* OVERLAPPING */
	    //  --pair--->
	    //<--this---
	    else if(pair.getMin() >= this.min)
		return 13;
	    // --pair--->
	    //  <---this--
	    else if(pair.getMin() < this.min)
		return 14;
	}
	return -1;
    }
    
    

    
    // -----this-----   ......if this dist is too FAR.....  -------other-----
    public boolean isTooFarFrom(DirectionalRange other){
	int dist = other.getMin() - this.max;
	if(dist > Constants.DM*10)
	    return true;
	return false;
    }
    
    
    public int getPaddingWhenTestingisWithinFacingFwd(DirectionalRange other, DirectionalRange thispair, DirectionalRange otherpair, boolean thispairIsInversionSegment){
	int padding = 0;
	if(thispairIsInversionSegment){
	    //int padding = 0;
	    int maxInversionLength = (other.getMax() > thispair.getMax() ? other.getMax() : thispair.getMax()) 
		- (other.getMin() < thispair.getMin() ? other.getMin() : thispair.getMin());
	    if(maxInversionLength < Constants.MIN_D)
		padding = Constants.MIN_D - maxInversionLength;
	}else{ // this DR is the inversion segment
	    int maxInversionLength = (otherpair.getMax() > this.getMax() ? otherpair.getMax() : this.getMax())
		- (otherpair.getMin() < this.getMin() ? otherpair.getMin() : this.getMin());
	    if(maxInversionLength < Constants.MIN_D)
		padding = Constants.MIN_D - maxInversionLength;
	}
	return padding;
    }



    public boolean isWithinFacingFwdINVERSION(DirectionalRange other, int padding){
	if(padding > Constants.MAX_CLUSTER_GAP){
	    //--this---> <--other--
	    if(this.isFwd() && !other.isFwd()){
		int dist = other.getMin() - this.max;
		if( dist > Constants.CLUSTER_OVERLAP 
		    && dist < padding
		    && this.min < other.getMin()
		    && this.max < other.getMax() )
		    return true;
	    }
	    return false;
	}else
	    return this.isWithinFacingFwd(other);
    }


    public boolean isWithinFacingFwdLenient(DirectionalRange other){
	if(Constants.MAX_CLUSTER_GAP < (Constants.MEDIAN/2)){
	    if(this.isFwd() && !other.isFwd()){
		int dist = other.getMin() - this.max;
		if( dist > Constants.CLUSTER_OVERLAP 
		    && dist < (Constants.MEDIAN/2)
		    && this.min < other.getMin()
		    && this.max < other.getMax() )
		    return true;
	    }
	    return false;
	}else
	    return this.isWithinFacingFwd(other);    
	
    }

    public boolean isWithinFacingFwd(DirectionalRange other){
	//--this---> <--other--
	if(this.isFwd() && !other.isFwd()){
	    int dist = other.getMin() - this.max;
	    if( dist > Constants.CLUSTER_OVERLAP 
		&& dist < Constants.MAX_CLUSTER_GAP  
		&& this.min < other.getMin()
		&& this.max < other.getMax() )
		return true;
	}
	return false;
    }


    // this is observed both in target (duplicative) insertion/transposition and inversion
    // trying to find DirectionalRange that are facing.
    //--DR1---> <--DR2--
    public boolean isWithinFacing(DirectionalRange other){
	//--this---> <--other--
	if(this.isWithinFacingFwd(other))
	    return true;
	/*if(this.isFwd() && !other.isFwd()){
	    int dist = other.getMin() - this.max;
	    if( dist > Constants.CLUSTER_OVERLAP 
		&& dist < Constants.MAX_CLUSTER_GAP  
		&& this.min < other.getMin()
		&& this.max < other.getMax() )
		return true;
		}*/
	//--other---> <--this--
	else if(other.isWithinFacingFwd(this))
	    return true;
	/*else if(!this.isFwd() && other.isFwd()){
	    int dist = this.min - other.getMax();
	    if( dist > Constants.CLUSTER_OVERLAP
		&& dist < Constants.MAX_CLUSTER_GAP
		&& other.getMin() < this.min
		&& other.getMax() < this.max )
		return true;
		}*/
	return false;
    }

    
    public boolean isBoundaryDefinerFwd(DirectionalRange other){
	//<----this----   ----other---->
	if(!this.isFwd() && other.isFwd()){
	    int dist = other.getMax() - this.min;
	    if( dist > Constants.MIN_TRANSPOSITION_TARGET 
		&& dist < Constants.MAX_TRANSPOSITION_TARGET )
		//&& this.min < other.getMin()     //removed these requirements to allow small transposition segment     |           |
		//&& this.max < other.getMax() )                                                                    --------other--->
		//                                                                                                        <-----this---------
		return true;
	}
	return false;
    }

    //trying to find DirectionalRange that are boundary defining of a target segment for transposition
    // <---DR1---    ----DR2--->
    public boolean isBoundaryDefiner(DirectionalRange other){
	//<----this----   ----other---->
	if(this.isBoundaryDefinerFwd(other))
	    return true;
	/*if(!this.isFwd() && other.isFwd()){
	    int dist = other.getMax() - this.min;
	    if( dist > Constants.MIN_TRANSPOSITION_TARGET 
		&& dist < Constants.MAX_TRANSPOSITION_TARGET
		&& this.min < other.getMin()
		&& this.max < other.getMax() )
		return true;
		}*/
	//<----other---   -----this---->
	else if(other.isBoundaryDefinerFwd(this))
	    return true;
	/*else if(this.isFwd() && !other.isFwd()){
	    int dist = this.max - other.getMin();
	    if( dist > Constants.MIN_TRANSPOSITION_TARGET
		&& dist < Constants.MAX_TRANSPOSITION_TARGET
		&& other.getMin() < this.min
		&& other.getMax() < this.max )
		return true;
		}*/
	return false;
    }

    public boolean isTandemDUPDirectional(DirectionalRange pair){
	int spanSize = 0;// = Constants.MIN_TANDUP_SPANNNIG_SIZE;
	if(!this.isFwd() && pair.isFwd()){
	    spanSize = pair.getMax() - this.min;
	
	    //if(spanSize >= Constants.MIN_TANDUP_SIZE
	    //   && spanSize <= Constants.MAX_TANDUP_SIZE)
	    if(spanSize <= Constants.MAX_TANDUP_SIZE)
		return true;
	}
	return false;
    }

    public boolean isSimpleDeletionDirectional(DirectionalRange pair){
	int dist = Constants.MIN_DELETION_SIZE;
	int spanSize = Constants.MIN_DELETION_SPANNING_SIZE_INCLUDING_CLUSTERS;
	if(this.isFwd() && !pair.isFwd()){
	    dist = pair.getMin() - this.max;
	    spanSize = pair.getMax() - this.min;
	}

	if(dist > Constants.MIN_DELETION_SIZE
	   && spanSize > Constants.MIN_DELETION_SPANNING_SIZE_INCLUDING_CLUSTERS)
	    return true;
	else 
	    return false;
    }
    
    public boolean isSimpleDeletionDirectional(DirectionalRange pair, boolean debug){
	
	int dist = Constants.MIN_DELETION_SIZE;
	int spanSize = Constants.MIN_DELETION_SPANNING_SIZE_INCLUDING_CLUSTERS;
	if(debug){
	    System.err.println("\t\t\t: Constants.MIN_DELETION_SPANNING :" + Constants.MIN_DELETION_SPANNING_SIZE_INCLUDING_CLUSTERS);
	    System.err.println("\t\t\t: Constants.MIN_DELETION_SIZE :" + Constants.MIN_DELETION_SIZE + "\tConstants.MIN_D:\t" + Constants.MIN_D );
	    System.err.println("\t\t\tPARAMETER: minDist = " + dist + "\tspanSize = " + spanSize);
	}
	if(this.isFwd() && !pair.isFwd()){
	    if(debug)
		System.err.println("\t\t\t------- ORIENTATION IS CORRECT!");
	    dist = pair.getMin() - this.max;
	    spanSize = pair.getMax() - this.min;
	    if(debug)
		System.err.println("\t\t\tCluser's dist = " + dist + "\tspanSize = " + spanSize);
	}else{
	    if(debug)
		System.err.println("\t\t\t------- ORIENTATION IS WRONG!!");
	}

	if(dist > Constants.MIN_DELETION_SIZE
	   && spanSize > Constants.MIN_DELETION_SPANNING_SIZE_INCLUDING_CLUSTERS){
	    if(debug)
		System.err.println("\t\t\t------- over min deletion SIZE. YEST!!");
	    return true;
	}else{ 
	    if(debug)
		System.err.println("\t\t\t------- SMALLER THAN MIN DELETION SIZE. NO!!");
	    return false;
	}
    }

    /* MUST BE CALLED from a pairing DirectionalRange! to be VALID*/
    // //----DR--->  -------some distance away   <------DR'(pair of DR)----
    // if the spanning gap is not repetative, then subsequent call to check coverage must be valided to confirm deletion.
    // if at least one of the DirectionalRange pair is from repetative Edge then this can also be a homologous recombination followed by duplicative transposition of the repeat.
    public boolean isSimpleDeletion(DirectionalRange pair){
	int dist = Constants.MIN_DELETION_SIZE;
	int spanSize = Constants.MIN_DELETION_SPANNING_SIZE_INCLUDING_CLUSTERS;
	
	return this.isSimpleDeletionDirectional(pair) || pair.isSimpleDeletionDirectional(this);
	
	//----this--->              <----pair----
	/*if(this.isFwd() && !pair.isFwd()){
	    dist = pair.getMin() - this.max;
	    spanSize = pair.getMax() - this.min;
	    }*/
	//----pair--->              <----this----
	/*else if(!this.isFwd() && pair.isFwd()){
	    dist = this.min - pair.getMax();
	    spanSize = this.max - pair.getMin();
	}

	if(dist > Constants.MIN_DELETION_SIZE
	   && spanSize > Constants.MIN_DELETION_SPANNING_SIZE_INCLUDING_CLUSTERS)
	    return true;
	    return false;*/
    }

    /*
     * @return: 
     * -1 : if not mergeable
     *  0 : NO changes in tip and stem --> "other" range is encompassed by "this"
     *  1 : if other is the tip --> "other" range is merged onto the tip side making itself the new updated tip
     *  2 : if other is the stem --> "other" range is merged onto the stem side making itself the new updated stem
     *  3 : if other is the tip & the stem --> "this" range is completely encompassed by "other" making itself" the new tip and stem.
     *  Tip is the interval that contains the end of the cluster(based on the direction)
     *  when HalfCluster spans 2 or more edges, 
     *
     *  ex) Interval     -----A------ ------B-------
     *      cluster             <-----------
     *  in this case the tip is interval A.
     *
     */
    public int mergeIfWithin(DirectionalRange other, int max_width){
	int nMin = this.min;
	int nMax = this.max;
	int tipStemVal = 0;
	if(this.isSameDirection(other)){
	    if(this.min > other.getMin()){
		nMin = other.getMin();
		if(!other.fwd)//[min tip, max stem]
		    tipStemVal = 1;
		else//[min stem, max tip]
		    tipStemVal = 2;
	    }
	    if(this.max < other.getMax()){
		nMax = other.getMax();
		if(other.fwd)//[max tip, min stem] : tipStemVal can only be either 0 (min not updated by "other" range) or 2. 
		    tipStemVal++;
		else//[max stem, min tip] : ]tipStemVal can only be either 0 (min not updated) or 1
		    tipStemVal = tipStemVal + 2;
	    }
	    if(nMax - nMin + 1 <= max_width){//if within max_width, update the value
		this.min = nMin;
		this.max = nMax;
		return tipStemVal;
	    }else{
		if(Constants.DEBUG)
		    System.err.println("[DirectionalRange.mergeIfWithin()]: mergedWidth(" +(nMax-nMin+1) +") exceeds max_width(" + max_width+ ")");
		return -1;
	    }
	}else
	    return -1;
    }
	

    
    public boolean isSameDirection(DirectionalRange other){
	if(this.fwd == other.fwd)
	    return true;
	return false;
    }

}


