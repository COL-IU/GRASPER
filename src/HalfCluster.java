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

//   HalfCluster keeps track of all the reads on one end of a full cluster. 
//   A full cluster always has two halfclusters for each end of the cluster. 
//   Reads belonging to the same end will be stored in one HalfCluster 
//   HalfCluster has reference to its pairing HalfCluster --> pair   
//   
//
//
//////////////////////////////////////////////////////////////////////
//   private int id;
//   private int curLinPos; --> this is used to store most recent 
//                              LinPos within activeClusterList. When 
//                              this is NOT in the activeClusterList, 
//                              it is set to its default value of -1
//
//   private HalfCluster pair;
//   private Read positionRead;
//
//   private ArrayList<Edge> edgeList;
//   private HashMap<Edge, EdgeCluster> edgeClusterHash; 
//
//   private ArrayList<DirectionalRange> directionalRangeList;
//   private int tipIndex; //index for accessing tip Edge in edgeList.
//   private int stemIndex; //index for accessing stem Edge in edgeList.
//////////////////////////////////////////////////////////////////////

public class HalfCluster{
    
    public HalfCluster(Read r1, Read r2){
	this(r1);
	this.setPair(new HalfCluster(r2));
	this.getPair().setPair(this);
    }
    
    //public HalfCluster(Read r, int nId){
    //	this(r);
    //	this.id = nextId;
    //	nextId++;
    //}

    public HalfCluster(Read r){
	
	this.checkRepetitiveDeletionAtTheEnd = false;

	this.directionalRangeList = null;
	this.tipIndex = 0;
	this.stemIndex = 0;

	this.id = nextId;
	nextId++;

	this.positionRead = r;
	this.edgeList = new ArrayList<Edge>();
	this.edgeList.add(r.getMappedEdge());
	this.edgeClusterHash = new HashMap<Edge, EdgeCluster>();
	//ArrayList<EdgeCluster> tmp = new ArrayList<EdgeCluster>();
	EdgeCluster ec = new EdgeCluster();
	ec.addRead(r);
	//tmp.add(ec);
	this.edgeClusterHash.put(r.getMappedEdge(), ec);
	this.pair = null;
	
	this.curLinPos = -1;
    }

    public void setCheckRepDeletionAgain(boolean check){
	this.checkRepetitiveDeletionAtTheEnd = check;
    }

    public boolean checkRepDeletionAgain(){
	return this.checkRepetitiveDeletionAtTheEnd;
    }

    public void setPair(HalfCluster other){
	this.pair = other;
    }

    public HalfCluster getPair(){
	return this.pair;
    }
    
    public ArrayList<DirectionalRange> getDirectionalRangeList(){
	return this.directionalRangeList;
    }

    public Read getPositionRead(){
	return this.positionRead;
    }



    //insertion sort
    public void sortEdgeListBasedOnMultiplicity(){
	for(int i = 1; i<this.edgeList.size(); i++){
	    Edge next = this.edgeList.get(i);
	    
	    int j = i;
	    while(j>0 && this.edgeList.get(j-1).getMultiplicity() > next.getMultiplicity()){
		this.edgeList.set(j,this.edgeList.get(j-1));
		j--;
	    }
	    this.edgeList.set(j,next);
	}
    }


    public void sortEdgeList(){
	int n = 0;//totalSize --> sum of #Intervals over edgeList

	ArrayList<Integer> listOfIndicies = new ArrayList<Integer>();
	HashMap<Integer, Edge> index2EdgeHash = new HashMap<Integer, Edge>();
	for(int i=0;i<this.edgeList.size();i++){
	    Edge e = this.edgeList.get(i);
	    n += e.getMultiplicity();
	    for(int j=0; j < e.getMultiplicity(); j++){
		int curThreadIndex = e.nthInterval(j).getThreadIndex();
		listOfIndicies.add(curThreadIndex);
		index2EdgeHash.put(new Integer(curThreadIndex),e);
	    }
	}
	Collections.sort(listOfIndicies);
	
    }
    
    public void removeDepletedRegionsThatMatchDeletionSign(ArrayList<DirectionalRange> deletionRangeList, int[] pairIndicies){
	DirectionalRange delDR = null;
	if(pairIndicies[0] == 0){
	    delDR = new DirectionalRange(this.directionalRangeList.get(pairIndicies[1]).getMax()
					 , this.pair.getDirectionalRangeList().get(pairIndicies[2]).getMin(), true, -1);
	}else if(pairIndicies[0] == 1){
	    delDR = new DirectionalRange(this.pair.getDirectionalRangeList().get(pairIndicies[1]).getMax()
					 , this.directionalRangeList.get(pairIndicies[2]).getMin(), true, -1);
	}else if(pairIndicies[0] == -1){
	    delDR = new DirectionalRange(this.directionalRangeList.get(pairIndicies[1]).getMax()
					 , this.pair.getDirectionalRangeList().get(pairIndicies[2]).getMin(), true, -1);
	}else if(pairIndicies[0] == -2){
	    delDR = new DirectionalRange(this.pair.getDirectionalRangeList().get(pairIndicies[1]).getMax()
					 , this.directionalRangeList.get(pairIndicies[2]).getMin(), true, -1);
	}
	
	for(int i=0;i<deletionRangeList.size();i++){
	    int containVal = delDR.contains(deletionRangeList.get(i));
	    if(containVal == 0){
		deletionRangeList.remove(i);
		i--;
	    }else if(containVal > 1){
		break;
	    }
	}
    }
    

    public void removeDepletedRegionsThatMatchDeletionSign(ArrayList<DirectionalRange> deletionRangeList, int deletionType){
		
	int nextJ = 0;
	for(int i=0; i<this.directionalRangeList.size();i++){
	    DirectionalRange curDR = this.directionalRangeList.get(i);
	    for(int j=nextJ; j<deletionRangeList.size();j++){
		int containVal = curDR.contains(deletionRangeList.get(j));
		if(containVal == 0){
		    deletionRangeList.remove(j);
		    j--;
		}else if(containVal > 1){
		    nextJ = j;
		    break;
		}
	    }
	}
    }

    
    /* SHOULDN"T BE USED. NOW A PART OF COMPUTEDIRECTIONALRANGES()*/
    public boolean sortEdgeList(Thread t){
	//size of the ArrayList should be --> number of Edges in edgeList
	//--->for each edge then it has an int[] size of  mutiplicity of that edge.
	//--------->Then each of the cell in int[] stores the threadIndex(intervalNumber)
	ArrayList<int[]> listOfIntervalIndicies = new ArrayList<int[]>();
	ECIndicies indicies = new ECIndicies(this.edgeList); //this stores bunch of indicies: current threadIndex, edgeIndex for locating an Edge from EdgeList, and EC indices.

	//load up listOfIntervalIndicies by running through Edges in edgeList.
	for(int i=0;i<this.edgeList.size();i++){
	    Edge e = this.edgeList.get(i);
	    int eMul = e.getMultiplicity();
	    listOfIntervalIndicies.add(new int[eMul]);
	    for(int j=0; j < eMul; j++){
		int curThreadIndex = e.nthInterval(j).getThreadIndex();
		listOfIntervalIndicies.get(i)[j] = curThreadIndex;
	    }
	}
	
	EdgeCluster eClu0 = null;
	EdgeCluster eClui = null;

	boolean sorted = false;
	boolean failed = false;//this gets turned to true if there is no path within Dm connecting Edges.
	while(true){
	    indicies.sort();//this sorts the intervals in an ascending order (i.e. 5'->3' direction of genome)
	    eClu0 = this.edgeClusterHash.get(indicies.getEdgeForIndex(this.edgeList, 0)); //since it's sorted, first interval (index 0) is the 5' most interval.
	    int minLinPos = eClu0.getLinPosMinForGivenEC(indicies.getECForEdgeForIndex(0)); //this is the current min linear position
	    boolean minLinPosDir = true;
	    if(minLinPos < 0){
		minLinPosDir = false;
		minLinPos = (minLinPos < 0 ? -1*minLinPos : minLinPos );
	    }
	    for(int i=1; i<indicies.length(); i++){
		//EdgeCluster eClu0 = this.edgeClusterHash.get(indicies.getEdgeForIndex(this.edgeList, 0));
		eClui = this.edgeClusterHash.get(indicies.getEdgeForIndex(this.edgeList, i));
		//if not within dm, we no longer need to look at more edges. Instead need to drop the lowest interval and feed the next.
		int curMax = 0;
		if((curMax = eClui.isRemoteWithinDmGivenMinPos(minLinPos, minLinPosDir, indicies.getECForEdgeForIndex(i), Constants.DM)) < 0){
		    if(!indicies.feedNext(this.edgeList))
			failed = true;
		    break;
		}else{//if within dm, 
		    //we check if it's the last edge
		    if(i == (indicies.length()-1))
			sorted = true;//then it's sorted
		    else//if not last edge
			;//we check next edge in the next iteration(edge)
		}
	    }
	    if(failed)
		break;
	    //IF SORTED: we change the order so that it reflects the correct ordering.
	    if(sorted){
		int[] edgeIndicies = indicies.getEdgeIndicies();
		ArrayList<Edge> tmp = new ArrayList<Edge>();
		for(int i=0; i<edgeIndicies.length; i++)
		    tmp.add(this.edgeList.get(edgeIndicies[i]));
		this.edgeList = tmp;
		break;
	    }
	}
	
	return sorted;
    }


    public ArrayList<Edge> getEdgeList(){
	return this.edgeList;
    }
    
    public HashMap<Edge, EdgeCluster> getEdgeClusterHash(){
	return this.edgeClusterHash;
    }

    public void updateEC(Edge e, int ec){
	this.edgeClusterHash.get(e).updateEC(ec);
    }

    public void setCurLinPos(int p){
	this.curLinPos = p;
    }
    
    public int getCurLinPos(){
	return this.curLinPos;
    }
    
    public void resetCurLinPos(){
	this.curLinPos = -1;
    }

    public boolean isWithinBoundary(int minPos){
	if(this.curLinPos >= minPos)
	    return true;
	return false;
    }

    //this is true IFF two halfClusters are within dm and their pairs as well
    public boolean isWithinDm(HalfCluster cl, int dm){
	//check (this , cl) and pair :(this.pair, cl.pair)
	if(this.isWithinDmSingle(cl, dm) && this.getPair().isWithinDmSingle(cl.getPair(), dm))
	    return true;
	if(Constants.DEBUG)
	    System.err.println("\tDO NOT COME WITHIN DM! FAIL");
	return false;
    }

    private boolean isWithinDmSingle(HalfCluster cl, int dm){
	if(Constants.DEBUG)
	    System.err.println("\t\t[isWithinDm] checking pair:\t" + this.getId() + "\t" + cl.getId());
	//checking distance between this and cl
	//in order check distance between this and cl, must check all edge-edge distance
	for(int i=0; i<this.edgeList.size(); i++){
	    for(int j=0; j<cl.getEdgeList().size(); j++){
		if( !edgeClusterHash.get(edgeList.get(i)).withinDm(cl.getEdgeClusterHash().get(cl.getEdgeList().get(j)), dm) ){
		    if(Constants.DEBUG)
			System.err.println("\t-->[isWithinDm] Fail here");
		    return false;
		}
	    }
	}
	if(Constants.DEBUG)
	    System.err.println("\t-->End within DM");
	return true;
    }


    private boolean isRemoteSameDirection(HalfCluster cl){
	boolean conflict = false;
	for(int i=0; i<this.edgeList.size();i++){
	    Edge e = this.edgeList.get(i);
	    EdgeCluster t_ec = this.edgeClusterHash.get(e);
	    EdgeCluster cl_ec = cl.getEdgeClusterHash().get(e);
	    if(cl_ec != null){
		if(cl_ec.minPosRead().isFwdAtEc(0) != t_ec.minPosRead().isFwdAtEc(0))
		    return false;
	    }else{
		for(int j=0;j<cl.getEdgeList().size();j++){
		    if(!e.areReadsSameDirectionClosestIntervals(cl.getEdgeList().get(j), t_ec.minPosRead(), cl.getEdgeClusterHash().get(cl.getEdgeList().get(j)).minPosRead()))
			return false;
		}
	    }
	}
	for(int j=0; j<cl.getEdgeList().size();j++){
	    Edge cl_e = cl.getEdgeList().get(j);
	    EdgeCluster cl_ec = cl.getEdgeClusterHash().get(cl_e);
	    EdgeCluster t_ec = this.edgeClusterHash.get(cl_e);
	    if(t_ec == null){
		for(int i=0;i<this.edgeList.size();i++){
		    if(!cl_e.areReadsSameDirectionClosestIntervals(this.edgeList.get(i), cl_ec.minPosRead(), this.edgeClusterHash.get(this.edgeList.get(i)).minPosRead()))
			return false;
		}
	    }
	}
	return true;
    }

    public void addAllFrom(HalfCluster cl){
	this.addAllFromSingle(cl);
	this.pair.addAllFromSingle(cl.getPair());
    }

    private void addAllFromSingle(HalfCluster cl){
	//need to check each edge;
	for(int i=0; i<cl.getEdgeList().size();i++){
	    EdgeCluster eClu = this.edgeClusterHash.get(cl.getEdgeList().get(i));
	    
	    //edge missing --> so need to add.
	    if(eClu == null){
		this.edgeList.add(cl.getEdgeList().get(i)); // add Edge
		this.edgeClusterHash.put(cl.getEdgeList().get(i), cl.getEdgeClusterHash().get(cl.getEdgeList().get(i))); //add the corresponding edgeCluster
	    }
	    //edge NOT missing --> need to add reads and update min/max
	    else{
		eClu.addReadsFromEdgeCluster(cl.getEdgeClusterHash().get(cl.getEdgeList().get(i)));
	    }
	}
    }
    
    //[10/8/14] : updated so we use the direction found in most recent visit.
    public boolean hasSameLeadingOrientationAs(HalfCluster other, int thisEC, int otherEC){
	if(Constants.DEBUG)
	    System.err.println("\tChecking leading orientation for:\t" + this.getId() + "\t" + other.getId());
	if( this.positionRead.isFwdAtEc(thisEC) == other.getPositionRead().isFwdAtEc(otherEC) ){//if( this.positionRead.isFwd() == other.getPositionRead().isFwd() )
	    if( this.pair.isRemoteSameDirection(other.getPair()) )
		return true;
	}
	if(Constants.DEBUG)
	    System.err.println("\tDO NOT HAVE SAME ORIENTATION! FAIL");
	return false;
    }

    public int getId(){
	return this.id;
    }

    public boolean isLeading(){
	return this.leading;
    }
    
    public void setLeadding(boolean l){
	this.leading = l;
	this.pair.setLeadingSingle(!l);
    }

    private void setLeadingSingle(boolean l){
	this.leading = l;
    }

    public boolean[] getInclusionArrayForSignificantEdges(){
	boolean[] inclusionArr = new boolean[this.edgeList.size()];
	int[] numReads = new int[this.edgeList.size()];
	int[] lengths = new int[this.edgeList.size()];
	int[] coverages = new int[this.edgeList.size()];
	
	for(int i=0; i<this.edgeList.size();i++){
	    EdgeCluster eclu = this.edgeClusterHash.get(this.edgeList.get(i));
	    lengths[i] = Math.abs(eclu.minPosRead().getGraphPos() - eclu.maxPosRead().getGraphPos()) + 1
		+ (Constants.DYNAMIC_READLEN ? (eclu.minPosRead().getLength()/2 + eclu.maxPosRead().getLength()/2) : Constants.READLEN );//Constants.READLEN;
	    numReads[i] = eclu.size();
	    coverages[i] = eclu.getTotalBases() / lengths[i];//numReads[i]*Constants.READLEN / lengths[i];
	    if(coverages[i] <= Constants.MIN_COVERAGE)
		inclusionArr[i] = false;
	    else
		inclusionArr[i] = true;
	}
	return inclusionArr;
    }
    

    /* returns all possible ranges */
    /*
    private ArrayList<DirectionalRange> computeDirectionalRanges(){
	ArrayList<DirectionalRange> list = new ArrayList<DirectionalRange>();
	EdgeCluster eclu = this.edgeClusterHash.get(this.edgeList.get(0));
	for(int i=0; i<this.edgeList.get(0).getMultiplicity();i++){
	    list.add( this.edgeList.get(0).getNthDirectionalRange(eclu.minPosRead(), eclu.maxPosRead(), i) );
	}
	
	if(this.edgeList.size() > 1){
	    for(int i=1;i<this.edgeList.size();i++){
		list = this.mergeNextEdgeToDirectionalRanges(list, i);
	    }
	}
	this.directionalRangeList = list;
	return list;
	}*/
    
    public void computeDirectionalRangesPair(){
	//System.err.println("\tPROCESSING\t1"  );
	this.computeDirectionalRanges();
	//System.err.println("\tPROCESSING\t2"  );
	this.pair.computeDirectionalRanges();
    }

    /* NEW VERSION returns the sorted list of DRs */
    private ArrayList<DirectionalRange> computeDirectionalRanges(){
	boolean debug = false;
	//if(this.id == 5357)
	//   debug = true;
	if(debug){
	    System.err.println("COMPUTE_DIRECTIONAL_RANGES for ClusterID:\t" + this.id);
	    Constants.printConstants();
	    System.err.println("========BEFORE COMPUTING");
	    this.printPlain();
	    System.err.println("========BEGIN COMPUTING");
	}
	ArrayList<DirectionalRange> list = new ArrayList<DirectionalRange>();
	
	//size of the ArrayList should be --> number of Edges in edgeList
	//--->for each edge then it has an int[] size of  mutiplicity of that edge.
	//--------->Then each of the cell in int[] stores the threadIndex(intervalNumber)
	ArrayList<int[]> listOfIntervalIndicies = new ArrayList<int[]>();
	ECIndicies indicies = new ECIndicies(this.edgeList); //this stores bunch of indicies: current threadIndex, edgeIndex for locating an Edge from EdgeList, and EC indices.

	/*if(debug)
	    System.err.println("ECIndicies length :\t" + indicies.length());
	*/
	//load up listOfIntervalIndicies by running through Edges in edgeList.
	for(int i=0;i<this.edgeList.size();i++){
	    Edge e = this.edgeList.get(i);
	    int eMul = e.getMultiplicity();
	    listOfIntervalIndicies.add(new int[eMul]);
	    for(int j=0; j < eMul; j++){
		int curThreadIndex = e.nthInterval(j).getThreadIndex();
		listOfIntervalIndicies.get(i)[j] = curThreadIndex;
	    }
	}
	
	EdgeCluster eClu0 = null;
	EdgeCluster eClui = null;

	boolean sorted = false;
	boolean failed = false;//this gets turned to true if there is no path within Dm connecting Edges.

	ArrayList<Edge> tmpEdgeList = new ArrayList<Edge>();

	while(!failed){
	    indicies.sort();//this sorts the intervals in an ascending order (i.e. 5'->3' direction of genome)
	    //if(debug)
	    //	System.err.println("After sorting:\t" + indicies.printCurIndicies());
	    
	    eClu0 = this.edgeClusterHash.get(indicies.getEdgeForIndex(this.edgeList, 0)); //since it's sorted, first interval (index 0) is the 5' most interval.
	    int minLinPos = eClu0.getLinPosMinForGivenEC(indicies.getECForEdgeForIndex(0)); //this is the current min linear position. this retusn + for fwd, - for reverse respect to genome at this interval.
	    boolean minLinPosDir = true;
	    if(minLinPos < 0){
		minLinPosDir = false;
		minLinPos = (minLinPos < 0 ? -1*minLinPos : minLinPos );
	    }
	    
	    if(indicies.length() == 1){
		int curMax = eClu0.getLinPosMaxForGivenEC(indicies.getECForEdgeForIndex(0));
		list.add(new DirectionalRange(minLinPos, curMax, minLinPosDir, 0, 0));
		//sorted = true;
		if(!indicies.feedNext(this.edgeList))
		    failed = true;
	    }else{
		//System.err.println("=======NEED TO SORT=====");
		for(int i=1; i<indicies.length(); i++){
		    //EdgeCluster eClu0 = this.edgeClusterHash.get(indicies.getEdgeForIndex(this.edgeList, 0));
		    eClui = this.edgeClusterHash.get(indicies.getEdgeForIndex(this.edgeList, i));
		    //if not within dm, we no longer need to look at more edges. Instead need to drop the lowest interval and feed the next.
		    int curMax = 0;
		    if((curMax = eClui.isRemoteWithinDmGivenMinPos(minLinPos, minLinPosDir, indicies.getECForEdgeForIndex(i), Constants.DM + Constants.READLEN)) < 0){
			if(!indicies.feedNext(this.edgeList)){
			    //System.err.println("=====>HERE");
			    failed = true;
			}
			break;
			
		    }else{//if within dm, 
			//we check if it's the last edge
			if(i == (indicies.length()-1)){
			    //System.err.println("=====>LAST EDGE");
			    //if first sorted, we save the order in tmp
			    if(!sorted){
				int[] edgeIndicies = indicies.getEdgeIndicies();
				for(int j=0; j<edgeIndicies.length; j++){
				    if(Constants.DEBUG)
					System.err.println("EdgeListSize:\t" + this.edgeList.size() + "\t" + "edgeIndicies["+j+"]="+ edgeIndicies[j]);
				    tmpEdgeList.add(this.edgeList.get(edgeIndicies[j]));
				}
			    }
			    sorted = true;//then it's sorted
			    int tipIndex = 0;
			    int stemIndex = indicies.length() - 1;
			    if(minLinPosDir){
				tipIndex = indicies.length() - 1;
				stemIndex = 0;
			    }
			    
			    list.add(new DirectionalRange(minLinPos, curMax, minLinPosDir, tipIndex, stemIndex));
			    
			    if(!indicies.feedAll(this.edgeList)){
				//System.err.println("=====>HERE2");
				failed = true;
				break;
			    }
			}
			else//if not last edge
			    ;//we check next edge in the next iteration(edge)
		    }
		}
		//System.err.println("====>HERE");
		
		if(failed){
		    if(sorted)
			this.edgeList = tmpEdgeList;//sets the ordering to the first sorted.
		    //System.err.println("======>HERE2");
		    break;
		}
	    }
	}
	//return sorted;
	this.directionalRangeList = list;
	return list;
    }

    




    /* this returns list of directional ranges that are satisfied by ALL edgeClusters within HalfClusters */
    private ArrayList<DirectionalRange> mergeNextEdgeToDirectionalRanges(ArrayList<DirectionalRange> list, int edgeIndex ){
	Edge e = this.edgeList.get(edgeIndex);
	ArrayList<DirectionalRange> nList = new ArrayList<DirectionalRange>();

	EdgeCluster eclu = this.edgeClusterHash.get(e);
	int nextJ = 0;

	for(int i=0;i<list.size();i++){
	    DirectionalRange tmpRange = list.get(i);
	    boolean merged = false;
	    for(int j=nextJ; j<e.getMultiplicity();j++){
		DirectionalRange tmpRange2 = e.getNthDirectionalRange(eclu.minPosRead(), eclu.maxPosRead(), j);
		int intervalThreadIndex = e.nthInterval(j).getThreadIndex();
		/* 
		 * mergedVal
		 * -1 : if not mergeable
		 *  0 : NO changes in tip and stem --> "other" range is encompassed by "this"
		 *  1 : if other is the tip --> "other" range is merged onto the tip side making itself the new updated tip
		 *  2 : if other is the stem --> "other" range is merged onto the stem side making itself the new updated stem
		 *  3 : if other is the tip & the stem --> "this" range is completely encompassed by "other" making itself" the new tip and stem.
		 *
		 * -1 : if not mergeable
		 *  1 : if other is the tip 
		 *  0 : if this is the tip
		 */
		int mergedVal = tmpRange.mergeIfWithin(tmpRange2, Constants.DM + Constants.READLEN); 
		//nextJ = j+1;
		if(mergedVal >= 0){
		    merged = true;
		    if(mergedVal == 1)
			tmpRange.setTipIntervalIndex(intervalThreadIndex);
		    else if(mergedVal == 2)
			tmpRange.setStemIntervalIndex(intervalThreadIndex);
		    else if(mergedVal == 3){
			tmpRange.setTipIntervalIndex(intervalThreadIndex);
			tmpRange.setStemIntervalIndex(intervalThreadIndex);
		    }
		    nList.add(tmpRange);
		    nextJ = j+1;
		    break;
		}
	    }
	    if(!merged){
		tmpRange.setTipIntervalIndex(-1);
		tmpRange.setStemIntervalIndex(-1);
	    }
	}
	return nList;
    }
    


    /*
     *  This returns the first linear position pairs and direction of the halfCluster in DirectionRange obj.
     */
    public DirectionalRange getDirectionalRange(){
	StringBuffer result = new StringBuffer();
	boolean[] inclusionArr = this.getInclusionArrayForSignificantEdges();
	int numEdges = 0;
	for(int i=0;i<inclusionArr.length;i++){
	    if(inclusionArr[i])
		numEdges++;
	}
	int first = 0;
	for(int i=0;i<this.edgeList.size();i++){
	    if(inclusionArr[i]){
		first = i;
		break;
	    }
	}
	EdgeCluster eclu = this.edgeClusterHash.get(this.edgeList.get(first));
	return this.edgeList.get(first).getFirstLinearPositionPairsAsDirectionalRange(eclu.minPosRead(), eclu.maxPosRead(), positionRead.isFwdAtEc(0));
    }

    public DirectionalRange getNthDirectionalRange(int n){
	return this.directionalRangeList.get(n);
    }

    public void printDirectionalRanges(){
	for(int i=0;i<this.directionalRangeList.size();i++){
	    System.err.println(i+ " :\t" + this.directionalRangeList.get(i).toString());
	}
    }
    
    /*
    public boolean isInversion(HalfCluster other){
	//index -->  i*pair.size + j
	//followerVals are defined as:
	// 0 :  ---this--->    ---pair--->
	// 1 :  <---this---    <---pair---
	// 2 :  <---pair---    <---this---
	// 3 :  ---pair--->    ---this---> 
	int[] followerValThis = new int[this.directionalRangeList.size()*this.pair.getDirectionalRangeList().size()];
	int[] followerValOther = new int[other.getDirectionalRangeList().size()*other.getPair().getDirectionalRangeList().size()];

	for(int i=0; i<this.directionalRangeList.size();i++){
	    for(int j=0; j<this.pair.getDirectionalRangeList().size();j++){
		followerValThis[i*this.pair.getDirectionalRangeList().size() + j] = this.getNthDirectionalRange(i).isFollower(this.pair.getNthDirectionalRange(j));
	    }
	}
	
	for(int i=0; i<other.getDirectionalRangeList().size();i++){
	    for(int j=0; j<other.getPair().getDirectionalRangeList().size();j++){
		followerValOther[i*other.getPair().getDirectionalRangeList().size() + j] = other.getNthDirectionalRange(i).isFollower(other.getPair().getNthDirectionalRange(j));
	    }
	}
	
	for(int i=0; i<followerValThis.length; i++){
	    for(int j=0; j<followerValOther.length; j++){
		if(followerValThis[i] == 0){
		    if(followerValOther[j] == 1){
			//checking for case   ----this----> <----other---              ----this.pair ---> <---other.pair----
			if(this.getNthDirectionalRange(i/this.pair.getDirectionalRangeList().size()).isWithinFacingFwd(other.getNthDirectionalRange(j/other.getPair().getDirectionalRangeList().size()))
			   && this.getNthDirectionalRange(i%this.pair.getDirectionalRangeList().size()).isWithinFacingFwd(other.getNthDirectionalRange(j%other.getPair().getDirectionalRangeList().size())))
			    {
				return  true;
			    }
		    }else if(followerValOther[j] == 2){
			//checking for case   ----this----> <---other.pair---          ----this.pair-----> <---other----
			if(this.getNthDirectionalRange(i/this.pair.getDirectionalRangeList().size()).isWithinFacingFwd(other.getNthDirectionalRange(j%other.getPair().getDirectionalRangeList().size()))
			   && this.getNthDirectionalRange(i%this.pair.getDirectionalRangeList().size()).isWithinFacingFwd(other.getNthDirectionalRange(j/other.getPair().getDirectionalRangeList().size())))
			   {
			       return  true;
			   }
		    }
		}else if(followerValThis[i] == 3){
		    if(followerValOther[j] == 1){
			//checking for case  ----this.pair---> <---other----           ----this----> <---other.pair---
			if(this.getNthDirectionalRange(i%this.pair.getDirectionalRangeList().size()).isWithinFacingFwd(other.getNthDirectionalRange(j/other.getPair().getDirectionalRangeList().size()))
			   && this.getNthDirectionalRange(i/this.pair.getDirectionalRangeList().size()).isWithinFacingFwd(other.getNthDirectionalRange(j%other.getPair().getDirectionalRangeList().size())))
			   {
			       return  true;
			   }
		    }else if(followerValOther[j] == 2){
			//checkinf for case ----this.pair---> <----other.pair---          ---this----> <---other---
			if(this.getNthDirectionalRange(i%this.pair.getDirectionalRangeList().size()).isWithinFacingFwd(other.getNthDirectionalRange(j%other.getPair().getDirectionalRangeList().size()))
			   && this.getNthDirectionalRange(i%this.pair.getDirectionalRangeList().size()).isWithinFacingFwd(other.getNthDirectionalRange(j%other.getPair().getDirectionalRangeList().size())))
			   {
			       return  true;
			   }
		    }
		}
	    }
	}
	return false;
    }

    */

    public ArrayList<DirectionalRange> getDirectionalRangeListForPairingType(int pairingType, boolean isOther){
	if(!isOther){//this
	    if(pairingType == 0 || pairingType == 1) //t
		return this.directionalRangeList;
	    else if(pairingType == 2 || pairingType == 3) //tp
		return this.getPair().getDirectionalRangeList();
	}else{//other
	    if(pairingType == 0 || pairingType == 2) //t
		return this.directionalRangeList;
	    else if(pairingType == 1 || pairingType == 3) //tp
		return this.getPair().getDirectionalRangeList();
	}
	System.err.println("INVALID pairingType: " + pairingType);
	System.exit(0);
	return null;
    }
    

    public ArrayList<DRWFCompareIndex> isWithinFacingForSpecificPairingType(HalfCluster other, int pairingType){
	return this.isWithinFacingForSpecificPairingType(other,pairingType, false);
    }
    //pairingType
    //0: T-O                                                                                                                                                              
    //1: T-OP                                                                                                                                                             
    //2: TP-O                                                                                                                                                             
    //3: TP-OP  
    //returns a position sorted list of DRWFIndicies                                                                                                                                                                                        
    //returns size zero list if none
    public ArrayList<DRWFCompareIndex> isWithinFacingForSpecificPairingType(HalfCluster other, int pairingType, boolean lenient){
	ArrayList<DRWFCompareIndex> drwfIndicies = new ArrayList<DRWFCompareIndex>();
	ArrayList<DirectionalRange> tdrl = this.getDirectionalRangeListForPairingType(pairingType, false);
	ArrayList<DirectionalRange> odrl = other.getDirectionalRangeListForPairingType(pairingType, true);
	
	for(int i=0; i<tdrl.size(); i++){
	    DirectionalRange tDR = tdrl.get(i);
	    for(int j=0; j<odrl.size(); j++){
		DirectionalRange oDR = odrl.get(j);
		int tDRFirstBit = -1;
		if(lenient){
		    if(tDR.isWithinFacingFwdLenient(oDR))
			tDRFirstBit = 1;
		    else if(oDR.isWithinFacingFwdLenient(tDR))
			tDRFirstBit = 0;
		    else
			tDRFirstBit = -1;
		}else{
		    if(tDR.isWithinFacingFwd(oDR))
			tDRFirstBit = 1;
		    else if(oDR.isWithinFacingFwd(tDR))
			tDRFirstBit = 0;
		    else
			tDRFirstBit = -1;
		}
		if(tDRFirstBit > -1)
		    drwfIndicies.add(new DRWFCompareIndex(tDR, i, oDR, j, (tDRFirstBit == 1 ? true : false), pairingType)); 
	    }
	}
	if(drwfIndicies.size() > 1)
	    Collections.sort(drwfIndicies, new DRWFPosComp());
	//if(drwfIndicies.size() == 0)
	//    return null;
	return drwfIndicies;
    }

    public ArrayList<DRWFCompareIndex> isWithinFacing(HalfCluster other){
	ArrayList<DRWFCompareIndex> drwfIndicies = new ArrayList<DRWFCompareIndex>();
	
	for(int i=0;i<this.directionalRangeList.size(); i++){
	    DirectionalRange tDR = this.getNthDirectionalRange(i);
	    
	    //checking tDR-oDR pairing
	    for(int j=0; j<other.getDirectionalRangeList().size(); j++){
		DirectionalRange oDR = other.getNthDirectionalRange(j);
		int tDRFirstBit = -1;
		if(tDR.isWithinFacingFwdLenient(oDR))
		    tDRFirstBit = 1;
		else if(oDR.isWithinFacingFwdLenient(tDR))
		    tDRFirstBit = 0;
		else
		    tDRFirstBit = -1;
		
		if(tDRFirstBit > -1)
		    drwfIndicies.add(new DRWFCompareIndex(tDR, i, oDR, j, (tDRFirstBit == 1 ? true : false), 0)); //0 for tDR-oDR pair
	    }
	    
	    //checking tDR-opDR pairing
	    for(int j=0; j<other.getPair().getDirectionalRangeList().size(); j++){
		DirectionalRange opDR = other.getPair().getNthDirectionalRange(j);
		int tDRFirstBit = -1;
		if(tDR.isWithinFacingFwdLenient(opDR))
		    tDRFirstBit = 1;
		else if(opDR.isWithinFacingFwdLenient(tDR))
		    tDRFirstBit = 0;
		else
		    tDRFirstBit = -1;

		if(tDRFirstBit > -1)
		    drwfIndicies.add(new DRWFCompareIndex(tDR, i, opDR, j, (tDRFirstBit == 1 ? true : false), 1));//1 for tDR-opDR pair
	    }
	}
	
	for(int i=0;i<this.getPair().getDirectionalRangeList().size(); i++){
	    DirectionalRange tpDR = this.getPair().getNthDirectionalRange(i);
	    
	    //checking tpDR-oDR pairing
	    for(int j=0; j<other.getDirectionalRangeList().size(); j++){
		DirectionalRange oDR = other.getNthDirectionalRange(j);
		int tpDRFirstBit = -1;
		if(tpDR.isWithinFacingFwdLenient(oDR))
		    tpDRFirstBit = 1;
		else if(oDR.isWithinFacingFwdLenient(tpDR))
		    tpDRFirstBit = 0;
		else
		    tpDRFirstBit = -1;

		if(tpDRFirstBit > -1)
		    drwfIndicies.add(new DRWFCompareIndex(tpDR, i, oDR, j, (tpDRFirstBit == 1 ? true : false), 2)); //2 for tpDR-oDR pair
	    }
	    
	    //checking tpDR-opDR pairing
	    for(int j=0; j<other.getPair().getDirectionalRangeList().size(); j++){
		DirectionalRange opDR = other.getPair().getNthDirectionalRange(j);
		int tpDRFirstBit = -1;
		if(tpDR.isWithinFacingFwdLenient(opDR))
		    tpDRFirstBit = 1;
		else if(opDR.isWithinFacingFwdLenient(tpDR))
		    tpDRFirstBit = 0;
		else
		    tpDRFirstBit = -1;
		
		if(tpDRFirstBit > -1)
		    drwfIndicies.add(new DRWFCompareIndex(tpDR, i, opDR, j, (tpDRFirstBit == 1 ? true : false), 3));//3 for tpDR-opDR pair
	    }
	}
	
	if(drwfIndicies.size() > 1)
	    Collections.sort(drwfIndicies, new DRWFPosComp());
	return drwfIndicies;
    }

    
    /*
    public ArrayList<int[]> isWithinFacing(HalfCluster other){
	ArrayList<int[]> listOfIndexPairForDirectionalRanges = new ArrayList<int[]> ();
	//directionalRangeList
	int nextJ = 0;
	for(int i=0;i<this.directionalRangeList.size();i++){
	    DirectionalRange curDR = this.getNthDirectionalRange(i);
	    for(int j=nextJ;j<other.getDirectionalRangeList().size();j++){
		if(curDR.isWithinFacing(other.getNthDirectionalRange(j))){
		    listOfIndexPairForDirectionalRanges.add(new int[] {i,j});
		    nextJ=j+1;
		    break;
		}else if(curDR.isTooFarFrom(other.getNthDirectionalRange(j))){
		    break;
		}
	    }
	}
	//if(listOfIndexPairForDirectionalRanges.size() > 0)
	return listOfIndexPairForDirectionalRanges;
	    //return null;
	    }*/


    //returns a sorted list of DRBDCompareIndex based on size and position
    //if none, returns zero-sized drbdIndicies
    public ArrayList<DRBDCompareIndex> isBoundaryDefiner(HalfCluster other){
	ArrayList<DRBDCompareIndex> drbdIndicies = new ArrayList<DRBDCompareIndex>();
	for(int i=0;i<this.directionalRangeList.size(); i++){
	    DirectionalRange tDR = this.getNthDirectionalRange(i);
	    
	    //checking tDR-oDR pairing
	    for(int j=0; j<other.getDirectionalRangeList().size(); j++){
		DirectionalRange oDR = other.getNthDirectionalRange(j);
		int tDRFirstBit = -1;
		if(tDR.isBoundaryDefinerFwd(oDR))
		    tDRFirstBit = 1;
		else if(oDR.isBoundaryDefinerFwd(tDR))
		    tDRFirstBit = 0;
		else
		    tDRFirstBit = -1;

		if(tDRFirstBit > -1)
		    drbdIndicies.add(new DRBDCompareIndex(tDR, i, oDR, j, (tDRFirstBit == 1 ? true : false), 0)); //0 for tDR-oDR pair
	    }
	    
	    //checking tDR-opDR pairing
	    for(int j=0; j<other.getPair().getDirectionalRangeList().size(); j++){
		DirectionalRange opDR = other.getPair().getNthDirectionalRange(j);
		int tDRFirstBit = -1;
		boolean tDRFirst = true;
		if(tDR.isBoundaryDefinerFwd(opDR)){
		    tDRFirstBit = 1;
		    tDRFirst = true;
		}else if(opDR.isBoundaryDefinerFwd(tDR)){
		    tDRFirstBit = 0;
		    tDRFirst = false;
		}else
		    tDRFirstBit = -1;

		if(tDRFirstBit > -1)
		    drbdIndicies.add(new DRBDCompareIndex(tDR, i, opDR, j, tDRFirst, 1));//1 for tDR-opDR pair
	    }
	}
	
	for(int i=0;i<this.getPair().getDirectionalRangeList().size(); i++){
	    DirectionalRange tpDR = this.getPair().getNthDirectionalRange(i);
	    
	    //checking tpDR-oDR pairing
	    for(int j=0; j<other.getDirectionalRangeList().size(); j++){
		DirectionalRange oDR = other.getNthDirectionalRange(j);
		int tpDRFirstBit = -1;
		if(tpDR.isBoundaryDefinerFwd(oDR))
		    tpDRFirstBit = 1;
		else if(oDR.isBoundaryDefinerFwd(tpDR))
		    tpDRFirstBit = 0;
		else
		    tpDRFirstBit = -1;
		if(tpDRFirstBit > -1)
		    drbdIndicies.add(new DRBDCompareIndex(tpDR, i, oDR, j, (tpDRFirstBit == 1 ? true : false), 2)); //2 for tpDR-oDR pair
	    }
	    
	    //checking tpDR-opDR pairing
	    for(int j=0; j<other.getPair().getDirectionalRangeList().size(); j++){
		DirectionalRange opDR = other.getPair().getNthDirectionalRange(j);
		int tpDRFirstBit = -1;
		boolean tpDRFirst = true;
		if(tpDR.isBoundaryDefinerFwd(opDR)){
		    tpDRFirstBit = 1;
		    tpDRFirst = true;
		}else if(opDR.isBoundaryDefinerFwd(tpDR)){
		    tpDRFirstBit = 0;
		    tpDRFirst = false;
		}else
		    tpDRFirstBit = -1;
		if(tpDRFirstBit > -1)
		    drbdIndicies.add(new DRBDCompareIndex(tpDR, i, opDR, j, tpDRFirst, 3));//3 for tpDR-opDR pair
	    }
	}
	
	//first sort based on size;
	if(drbdIndicies.size() > 1){
	    Collections.sort(drbdIndicies, new SizeComp());
	
	    //second we sort subList in MAX_D size range
	    int start = 0;
	    if(drbdIndicies.size() > 0){
		int uptoSize = drbdIndicies.get(0).getSize() + Constants.MAX_D;
		
		for(int i=1;i<drbdIndicies.size();i++){
		    if(drbdIndicies.get(i).getSize() >= uptoSize){
			uptoSize = drbdIndicies.get(i).getSize() + Constants.MAX_D;
			Collections.sort(drbdIndicies.subList(start,i), new PositionComp());
			start = i;
		    }else if( i == (drbdIndicies.size()-1) )
			Collections.sort(drbdIndicies.subList(start,i), new PositionComp());
		}
	    }
	}
	
	return drbdIndicies;
    }

    public String getEventType(HalfCluster other, boolean debug){
	
	System.err.println("****" + this.getId() + "\t :\t" + other.getId());
	
	System.err.println("---- T DRs");
	this.printDirectionalRanges();
	
	System.err.println("---- TP DRs");
	this.pair.printDirectionalRanges();
	
	System.err.println("---- O DRs");
	other.printDirectionalRanges();
	
	System.err.println("---- OP DRs");
	other.getPair().printDirectionalRanges();
	
	ArrayList<DRBDCompareIndex> drbdIndicies = this.isBoundaryDefiner(other);

	//key: pairingType, Value: sorted DRWFComapreIndex by genomic position.
	HashMap<Integer, ArrayList<DRWFCompareIndex>> drwfIndexHash = new HashMap<Integer, ArrayList<DRWFCompareIndex>>();
	
	DirectionalRange tDR = null;
	DirectionalRange tpDR = null;
	DirectionalRange oDR = null;
	DirectionalRange opDR = null;
	
	for(int i=0; i<drbdIndicies.size(); i++){
	    DRBDCompareIndex drbdIndex = drbdIndicies.get(i);
	    //0: T-O
	    //1: T-OP
	    //2: TP-O
	    //3: TP-OP  
	    int drbdPT = drbdIndex.pairingType();
	    int drwfPT = drbdIndex.getInvertedPT();
	    System.err.println("drbdPT:drwfPT\t" + drbdPT + ":" + drwfPT );
	    
	    if(drwfIndexHash.get(new Integer(drwfPT)) == null){
		drwfIndexHash.put(new Integer(drwfPT), this.isWithinFacingForSpecificPairingType(other, drwfPT, false));
	    }
	    ArrayList<DRWFCompareIndex> drwfIndicies = drwfIndexHash.get(new Integer(drwfPT));
	    
	    if(drwfIndicies.size() > 0){
		DRWFCompareIndex drwfIndex = drwfIndicies.get(0); // gets the first one since it's sorted
		//tDR = this.getNthDirectionalRange(drbdIndex.tDRIndex());
		//oDR = other.getNthDirectionalRange(drbdIndex.oDRIndex());
		//tpDR = this.pair.getNthDirectionalRange(drwfIndex.tDRIndex());
		//opDR = other.getPair().getNthDirectionalRange(drwfIndex.oDRIndex());
		
		if(drbdPT == 0){ //drwfPT is 3
		    System.err.println("0-3 HERE\t" + drbdPT + ":" + drwfPT);
		    tDR = this.getNthDirectionalRange(drbdIndex.tDRIndex());
		    oDR = other.getNthDirectionalRange(drbdIndex.oDRIndex());
		    tpDR = this.pair.getNthDirectionalRange(drwfIndex.tDRIndex());
		    opDR = other.getPair().getNthDirectionalRange(drwfIndex.oDRIndex());
		    System.err.println("tDR:\t" + tDR);
		    System.err.println("oDR:\t" + oDR);
		    System.err.println("tpDR:\t" + tpDR);
		    System.err.println("opDR:\t" + opDR);
		    // <--this---   ----other-->
		    // -------IS-ELEMENT--------
		    if(drbdIndex.tDRFirst()){
			// ---this.p--> <---other.p---
			// ------------|-------------
			// insertion site "\t|\t"
			if(drwfIndex.tDRFirst()){
			    System.out.println("TRANSPOSTION\t"+"1-1"+"\t" + this.getNumReadsDRSize() +other.getNumReadsDRSize()
					       + tDR.toString() +"\t"+ oDR.toString() + "\t|\t" + tpDR.toString() +"\t"+ opDR.toString()
					       + "\t" + this.pair.getNumReadsDRSize() + other.getPair().getNumReadsDRSize());
			    return "1-1";//ASSIGN_TRANSPOSITION_1-1;
			}
			// ---other.p--> <---this.p---
			// -------------|-------------
			// insertion site "\t|\t"
			else{
			    
			    System.out.println("TRANSPOSTION(I)\t"+"5-1"+"\t" + this.getNumReadsDRSize() + other.getNumReadsDRSize()
					       + tDR.toString() +"\t" + oDR.toString() + "\t|\t" + opDR.toString()+"\t" + tpDR.toString()
					       + "\t" + other.getPair().getNumReadsDRSize() + this.pair.getNumReadsDRSize());
			    return "5-1";//ASSIGN_TRANSPOSITION_INVERSION_5-1;
			}
			
		    }
		    // <--other---   ----this-->
		    // -------IS-ELEMENT--------
		    else{
			// ---this.p--> <---other.p---
			// ------------|-------------
			// insertion site "\t|\t"
			if(drwfIndex.tDRFirst()){
			    System.out.println("TRANSPOSTION(I)\t"+"6-1"+"\t" + other.getNumReadsDRSize() + this.getNumReadsDRSize()
					       + oDR.toString() +"\t"+ tDR.toString() + "\t|\t" + tpDR.toString()+"\t" + opDR.toString()
					       + "\t" + this.pair.getNumReadsDRSize() + other.getPair().getNumReadsDRSize());
			    return "6-1";//ASSIGN_TRANSPOSITION_INVERTED_6-1;
			}
			// ---other.p--> <---this.p---
			// -------------|-------------
			// insertion site "\t|\t"
			else{
			    System.out.println("TRANSPOSTION\t"+"3-1"+"\t" + other.getNumReadsDRSize() + this.getNumReadsDRSize() 
					       + oDR.toString() +"\t"+ tDR.toString() + "\t|\t" + opDR.toString()+"\t" + tpDR.toString()
					       + "\t" + other.getPair().getNumReadsDRSize() + this.pair.getNumReadsDRSize());
			    return "3-1";//ASSIGN_TRANSPOSITION_3-1;
			}
		    }
		}
		else if(drbdPT == 1){//drwfPT is 2
		    tDR = this.getNthDirectionalRange(drbdIndex.tDRIndex());
		    opDR = other.getPair().getNthDirectionalRange(drbdIndex.oDRIndex());
		    tpDR = this.pair.getNthDirectionalRange(drwfIndex.tDRIndex());
		    oDR = other.getNthDirectionalRange(drwfIndex.oDRIndex());
		    // <--this---   ----other.p-->
		    // -------IS-ELEMENT--------
		    if(drbdIndex.tDRFirst()){
			
			// ---this.p--> <---other----
			// ------------|-------------
			// insertion site "\t|\t"
			if(drwfIndex.tDRFirst()){
			    System.out.println("TRANSPOSTION\t"+"2-1"+"\t" + this.getNumReadsDRSize() + other.getPair().getNumReadsDRSize()
					       + tDR.toString() +"\t"+ opDR.toString() + "\t|\t" + tpDR.toString() +"\t"+ oDR.toString()
					       + "\t" + this.pair.getNumReadsDRSize() + other.getNumReadsDRSize());
			    return "2-1";//ASSIGN_TRANSPOSITION_2-1;
			}
			// ---other--> <---this.p---
			// -----------|-------------
			// insertion site "\t|\t"
			else{
			    System.out.println("TRANSPOSTION(I)\t"+"5-2"+"\t" + this.getNumReadsDRSize() + other.getPair().getNumReadsDRSize() 
					       + tDR.toString() + "\t" +opDR.toString() + "\t|\t" + oDR.toString() +"\t"+ tpDR.toString()
					       + "\t" + other.getNumReadsDRSize() + this.pair.getNumReadsDRSize());
			    return "5-2";//ASSIGN_TRANSPOSITION_INVERSION_5-2;
			}
		    }
		    // <--other.p--   ---this-->
		    // -------IS-ELEMENT--------
		    else{
			
			// ---this.p--> <---other---
			// ------------|-------------
			// insertion site "\t|\t"
			if(drwfIndex.tDRFirst()){
			    System.out.println("TRANSPOSTION(I)\t"+"6-3"+"\t" +other.getPair().getNumReadsDRSize() + this.getNumReadsDRSize() 
					       + opDR.toString() +"\t"+ tDR.toString() + "\t|\t" + tpDR.toString() +"\t"+ oDR.toString()
					       + "\t" + this.pair.getNumReadsDRSize() + other.getNumReadsDRSize());
			    return "6-3";//ASSIGN_TRANSPOSITION_INVERTED_6-3;
			}
			// ---other--> <---this.p---
			// -----------|-------------
			// insertion site "\t|\t"
			else{
			    System.out.println("TRANSPOSTION\t"+"4-2"+"\t" + other.getPair().getNumReadsDRSize() + this.getNumReadsDRSize()
					       + opDR.toString() +"\t"+ tDR.toString() + "\t|\t" + oDR.toString() +"\t"+ tpDR.toString()
					       + "\t" + other.getNumReadsDRSize() + this.pair.getNumReadsDRSize());
			    return "4-2";//ASSIGN_TRANSPOSITION_4-2;
			}
		    }
		}
		else if(drbdPT == 2){//drwfPT is 1
		    tpDR = this.pair.getNthDirectionalRange(drbdIndex.tDRIndex());
		    oDR = other.getNthDirectionalRange(drbdIndex.oDRIndex());
		    tDR = this.getNthDirectionalRange(drwfIndex.tDRIndex());
		    opDR = other.getPair().getNthDirectionalRange(drwfIndex.oDRIndex());
		    // <--this.p---   ----other-->
		    // -------IS-ELEMENT--------
		    if(drbdIndex.tDRFirst()){
			
			// ---this----> <---other.p---
			// ------------|-------------
			// insertion site "\t|\t"
			if(drwfIndex.tDRFirst()){
			    System.out.println("TRANSPOSTION\t"+"4-1"+"\t" + this.pair.getNumReadsDRSize() + other.getNumReadsDRSize()
					       + tpDR.toString() +"\t"+ oDR.toString() + "\t|\t" + tDR.toString() +"\t"+ opDR.toString()
					       + "\t" + this.getNumReadsDRSize() + other.getPair().getNumReadsDRSize());
			    return "4-1";//ASSIGN_TRANSPOSITION_4-1;
			}
			// ---other.p--> <---this---
			// -------------|-------------
			// insertion site "\t|\t"
			else{
			    System.out.println("TRANSPOSTION(I)\t"+"5-3"+"\t" + this.pair.getNumReadsDRSize() + other.getNumReadsDRSize()
					       + tpDR.toString() + "\t" +oDR.toString() + "\t|\t" + opDR.toString() +"\t"+ tDR.toString()
					       + "\t" + other.getPair().getNumReadsDRSize() + this.getNumReadsDRSize());
			    return "5-3";//ASSIGN_TRANSPOSITION_INVERSION_5-3;
			}
		    }
		    // <--other---   ----this.p-->
		    // -------IS-ELEMENT--------
		    else{
			// ---this--> <---other.p---
			// ----------|-------------
			// insertion site "\t|\t"
			if(drwfIndex.tDRFirst()){
			    System.out.println("TRANSPOSTION(I)\t"+"6-2"+"\t" +other.getNumReadsDRSize() +  this.pair.getNumReadsDRSize()
					       + oDR.toString()+"\t" + tpDR.toString() + "\t|\t" + tDR.toString() +"\t"+ opDR.toString()
					       + "\t" + this.getNumReadsDRSize()  + other.getPair().getNumReadsDRSize());
			    return "6-2";//ASSIGN_TRANSPOSITION_INVERTED_6-2;
			}
			// ---other.p--> <---this----
			// ------------|-------------
			// insertion site "\t|\t"
			else{
			    System.out.println("TRANSPOSTION\t"+"2-2"+"\t" + other.getNumReadsDRSize() + this.pair.getNumReadsDRSize() 
					       + oDR.toString() +"\t"+ tpDR.toString() + "\t|\t" + opDR.toString() +"\t"+ tDR.toString()
					       + "\t" + other.getPair().getNumReadsDRSize() + this.getNumReadsDRSize());
			    return "2-2";//ASSIGN_TRANSPOSITION_2-2;
			}
		    }
		}
		else if(drbdPT == 3){//drwfPT is 0
		    tpDR = this.pair.getNthDirectionalRange(drbdIndex.tDRIndex());
		    opDR = other.getPair().getNthDirectionalRange(drbdIndex.oDRIndex());
		    tDR = this.getNthDirectionalRange(drwfIndex.tDRIndex());
		    oDR = other.getNthDirectionalRange(drwfIndex.oDRIndex());
		    // <--this.p---   ----other.p-->
		    // ----------IS-ELEMENT---------
		    if(drbdIndex.tDRFirst()){
			// ---this---> <---other---
			// -----------|-------------
			// insertion site "\t|\t"
			if(drwfIndex.tDRFirst()){
			    System.out.println("TRANSPOSTION\t"+"3-2"+"\t" + this.pair.getNumReadsDRSize() + other.getPair().getNumReadsDRSize() 
					       + tpDR.toString() +"\t"+ opDR.toString() + "\t|\t" + tDR.toString() +"\t"+ oDR.toString() 
					       + "\t" + this.getNumReadsDRSize() +other.getNumReadsDRSize());
			    return "3-2";//ASSIGN_TRANSPOSITION_3-2;
			}
			// ---other---> <---this---
			// ------------|-------------
			// insertion site "\t|\t"
			else{
			    System.out.println("TRANSPOSTION(I)\t"+"5-4"+"\t" + this.pair.getNumReadsDRSize() + other.getPair().getNumReadsDRSize()
					       + tpDR.toString() +"\t"+ opDR.toString() + "\t|\t" + oDR.toString() +"\t"+ tDR.toString()
					       + "\t" + other.getNumReadsDRSize() + this.getNumReadsDRSize() );
			    return "5-4";//ASSIGN_TRANSPOSITION_INVERSION_5-4;
			}
		    }
		    // <--other.p---   ----this.p-->
		    // ---------IS-ELEMENT----------
		    else{
			// ---this--> <---other---
			// ----------|------------
			// insertion site "\t|\t"
			if(drwfIndex.tDRFirst()){
			    System.out.println("TRANSPOSTION(I)\t"+"6-4"+"\t" + other.getPair().getNumReadsDRSize() + this.pair.getNumReadsDRSize()
					       + opDR.toString() +"\t"+ tpDR.toString() + "\t|\t" + tDR.toString() +"\t"+ oDR.toString()
					       + "\t" + this.getNumReadsDRSize() + other.getNumReadsDRSize());
			    return "6-4";//ASSIGN_TRANSPOSITION_INVERTED_6-4;
			}
			// ---other---> <---this---
			// ------------|-------------
			// insertion site "\t|\t"
			else{
			    System.out.println("TRANSPOSTION\t"+"1-2"+"\t" + other.getPair().getNumReadsDRSize() + this.pair.getNumReadsDRSize()
					       + opDR.toString() +"\t"+ tpDR.toString() + "\t|\t" + oDR.toString() +"\t"+ tDR.toString()
					       + "\t" + other.getNumReadsDRSize() + this.getNumReadsDRSize());
			    return "1-2";//ASSIGN_TRANSPOSITION_1-2;
			}
		    }
		}
	    }
	}
	

	
	/*ArrayList<DRWFCompareIndex> drwfIndiciesZero = drwfIndexHash.get(new Integer(0));
	  ArrayList<DRWFCompareIndex> drwfIndiciesThree = drwfIndexHash.get(new Integer(3));
	  
	  if(drwfIndiciesZero == null){
	  drwfIndexHash.put(new Integer(0), this.isWithinFacingForSpecificPairing(other, 0));
	  drwfIndiciesZero = drwfIndexHash.get(new Integer(0));
	  }
	*/
	
	int curInversionMinSize = Constants.GENOMELEN;
	String curInversionString = null;
	String inversionReturnVal = null;
	//0-3 pair
	ArrayList<DRWFCompareIndex> drwfIndiciesZero = this.isWithinFacingForSpecificPairingType(other, 0);
	ArrayList<DRWFCompareIndex> drwfIndiciesThree = null;
	if(drwfIndiciesZero.size() > 0){
	    
	    drwfIndiciesThree = this.isWithinFacingForSpecificPairingType(other, 3);
	    
	    if(drwfIndiciesThree.size() > 0){
		
		for(int i=0;i<drwfIndiciesZero.size(); i++){
		    for(int j=0; j<drwfIndiciesThree.size(); j++){
			DirectionalRange[] drs = drwfIndiciesZero.get(i).isInvertionPair(drwfIndiciesThree.get(j), this, other);
			if(drs != null){
			    int inversionType = drs[4].getInversionType();
			    int iSize = drs[3].getMin() - drs[0].getMax();
			    if(iSize < curInversionMinSize){
				curInversionString = "INVERSION\tI-" + inversionType + "\t" + this.getPartInversionString(inversionType, other, true)
				    + drs[0].toString() + "\t" + drs[1].toString() + "\t|\t"+ drs[2].toString() + "\t" + drs[3].toString()
				    + "\t" + this.getPartInversionString(inversionType, other, false);
				inversionReturnVal = "I-" + inversionType;
				curInversionMinSize = iSize;
			    }		       
			}
		    }
		}
	    }
	}
	
	//1-2 pair
	ArrayList<DRWFCompareIndex> drwfIndiciesOne = this.isWithinFacingForSpecificPairingType(other, 1);
	ArrayList<DRWFCompareIndex> drwfIndiciesTwo = null;
	
	if(drwfIndiciesOne.size() > 0){
	 
	    drwfIndiciesTwo = this.isWithinFacingForSpecificPairingType(other, 2);
	    
	    if(drwfIndiciesTwo.size() > 0){
		
		for(int i=0;i<drwfIndiciesOne.size(); i++){
		    for(int j=0; j<drwfIndiciesTwo.size(); j++){
			DirectionalRange[] drs = drwfIndiciesOne.get(i).isInvertionPair(drwfIndiciesTwo.get(j), this, other);
			if(drs != null){
			    int inversionType = drs[4].getInversionType();
			    int iSize =  drs[3].getMin() - drs[0].getMax();
			    if(iSize < curInversionMinSize){
				curInversionString = "INVERSION\tI-" + inversionType + "\t" + this.getPartInversionString(inversionType, other, true)
				    + drs[0].toString() + "\t" + drs[1].toString() + "\t|\t"+ drs[2].toString() + "\t" + drs[3].toString()
				    + "\t" + this.getPartInversionString(inversionType, other, false);
				inversionReturnVal = "I-" + inversionType;
				curInversionMinSize = iSize;
			    }		       
			}
		    }
		}
	    }
	}
	
	if(curInversionString != null){
	    System.out.println(curInversionString);
	    return inversionReturnVal;
	}
	return null;
	
    }

    private String getPartInversionString(int inversion_type, HalfCluster other, boolean first){
	int it = inversion_type;
	if(it == 1){
	    if(first)
		return this.getNumReadsDRSize() + other.getNumReadsDRSize();
	    else
		return this.getPair().getNumReadsDRSize() + other.getPair().getNumReadsDRSize();
	}else if(it == 2){
	    if(first)
		return this.getPair().getNumReadsDRSize() + other.getPair().getNumReadsDRSize();
	    else
		return this.getNumReadsDRSize() + other.getNumReadsDRSize();
	}else if(it == 3){
	    if(first)
		return other.getNumReadsDRSize() + this.getNumReadsDRSize();
	    else
		return other.getPair().getNumReadsDRSize() + this.getPair().getNumReadsDRSize();
	}else if(it == 4){
	    if(first)
		return other.getPair().getNumReadsDRSize() + this.getPair().getNumReadsDRSize();
	    else
		return other.getNumReadsDRSize() + this.getNumReadsDRSize();
	}else if(it == 5){
	    if(first)
		return this.getNumReadsDRSize() + other.getPair().getNumReadsDRSize();
	    else
		return this.getPair().getNumReadsDRSize() + other.getNumReadsDRSize();
	}else if(it == 6){
	    if(first)
		return this.getPair().getNumReadsDRSize() + other.getNumReadsDRSize();
	    else
		return this.getNumReadsDRSize() + other.getPair().getNumReadsDRSize();
	}else if(it == 7){
	    if(first)
		return other.getPair().getNumReadsDRSize() + this.getNumReadsDRSize();
	    else
		return other.getNumReadsDRSize() + this.getPair().getNumReadsDRSize();
	}else if(it == 8){
	    if(first)
		return other.getNumReadsDRSize() + this.getPair().getNumReadsDRSize();
	    else
		return other.getPair().getNumReadsDRSize() + this.getNumReadsDRSize();
	}
	return null;
    }



    /*
    public String getEventString(HalfCluster other, String eventType){
	int typeID = ClusterApp.getEventTypeID(eventType);
	String eventStr;
	StringBuffer tmp = new StringBuffer();
	if(typeID == 0 || typeID == 1){
	    tmp.append("TRANSPOSITION\t");
	    if(eventType.equals("1-1"))
		tmp.append(this.)
	}else if(typeID = 2){
	    tmp.append("INVERSION\t");
	}
	    
	}*/

    
    public String getNumReadsDRSize(){
	return this.numReads() + "\t" + this.directionalRangeList.size() + "\t";
    }


    public String getTandemDuplicationEventString(int[] pairIndicies){
	String tdStr;
	if(pairIndicies[0] == 0){
	    tdStr = "TANDEM_DUP\tT-1\t" + this.getNumReadsDRSize() + this.pair.getNumReadsDRSize() 
		+ this.directionalRangeList.get(pairIndicies[1]).toString() 
		+ "\t" + this.pair.getDirectionalRangeList().get(pairIndicies[2]).toString();
	}else if(pairIndicies[0] ==1){
	    tdStr = "TANDEM_DUP\tT-2\t" + this.pair.getNumReadsDRSize() + this.getNumReadsDRSize()
		+ this.pair.getDirectionalRangeList().get(pairIndicies[1]).toString() 
		+ "\t" + this.directionalRangeList.get(pairIndicies[2]).toString();
	}else
	    tdStr = "ERROR!!!!!!!!!!";
	
	System.out.println(tdStr);
	return tdStr;
    }

    public String getDeletionEventString(int[] pairIndicies){
	String delStr;
	if(pairIndicies[0] == 0){
	    delStr = "DELETION\tD-1\t" + this.getNumReadsDRSize() + this.pair.getNumReadsDRSize() 
		+ this.directionalRangeList.get(pairIndicies[1]).toString() 
		+ "\t" + this.pair.getDirectionalRangeList().get(pairIndicies[2]).toString();
	}else if(pairIndicies[0] == 1){
	    delStr = "DELETION\tD-2\t" + this.pair.getNumReadsDRSize() + this.getNumReadsDRSize()
		+ this.pair.getDirectionalRangeList().get(pairIndicies[1]).toString() 
		+ "\t" + this.directionalRangeList.get(pairIndicies[2]).toString();
	}else if(pairIndicies[0] == -1){
	    delStr = "DELETION\tR-1\t" + this.getNumReadsDRSize()  + this.pair.getNumReadsDRSize()
		+ this.directionalRangeList.get(pairIndicies[1]).toString() 
		+ "\t" + this.pair.getDirectionalRangeList().get(pairIndicies[2]).toString();
	}else if(pairIndicies[0] == -2){
	    delStr = "DELETION\tR-2\t" + this.pair.getNumReadsDRSize() + this.getNumReadsDRSize()
		+ this.pair.getDirectionalRangeList().get(pairIndicies[1]).toString() 
		+ "\t" + this.directionalRangeList.get(pairIndicies[2]).toString();
	}else
	    delStr = "ERROR!!!!!!!!!!!";

	System.out.println(delStr);
	return delStr;
    }


    //return array size of 3. [0]: 0 if this-pair deletion, 1 if pair-this deletion, -1 this-pair-rep deletion, -2 pair-this-rep deletion
    public int[] isSimpleDeletion(Thread t){
	//ArrayList<int[]> listOfIndexPairForDirectionalRanges = new ArrayList<int[]>();
	
	boolean debug = false;
	/*if(this.id == 8493 || this.id == 8494){
	   debug = true;
	   System.out.println("============14631||14632 FOUND==== CHECKING IF DELETION");
	   }*/
	boolean foundDeletion_this_pair = false;
	boolean foundDeletion_pair_this = false;

	boolean foundDeletion_this_pair_rep = false;
	boolean foundDeletion_pair_this_rep = false;

	int[] ij_this_pair = new int[2];
	int[] ij_pair_this = new int[2];

	int deletionSize_this_pair = 0;
	int deletionSize_pair_this = 0;

	int deletionSize_rep = Constants.GENOMELEN;

	for(int i=0; i<this.directionalRangeList.size();i++){
	    DirectionalRange rangeI = this.directionalRangeList.get(i);
	    if(rangeI.isFwd()){
		for(int j=0;j<this.pair.getDirectionalRangeList().size();j++){
		    DirectionalRange rangeJ = this.pair.getDirectionalRangeList().get(j);
		    if(!rangeJ.isFwd()){
			int tmpSize = rangeJ.getMin() - rangeI.getMax() - 1;
			boolean update = true;
			if(foundDeletion_this_pair){
			    if(tmpSize <= deletionSize_this_pair)
				update = false;
			}
			if(update){
			    if( (rangeI.getMax() < rangeJ.getMin())
				&& rangeI.isSimpleDeletionDirectional(rangeJ) ){
				int lackingCoverageVal = t.isPathLackingCoverage(new IntervalPos(t, rangeI.getMax()), 
										 new IntervalPos(t, rangeJ.getMin()));
				if(lackingCoverageVal == 1){
				    foundDeletion_this_pair = true;
				    deletionSize_this_pair = tmpSize;
				    ij_this_pair[0] = i;
				    ij_this_pair[1] = j;
				    if(foundDeletion_this_pair_rep){
					this.checkRepetitiveDeletionAtTheEnd = false;
					this.pair.setCheckRepDeletionAgain(false);
					foundDeletion_this_pair_rep = false;
				    }
				}else if(!foundDeletion_this_pair && lackingCoverageVal == -1){
				    if(tmpSize < deletionSize_rep){
					this.checkRepetitiveDeletionAtTheEnd = true;
					foundDeletion_this_pair_rep = true;
					this.pair.setCheckRepDeletionAgain(true);
					ij_this_pair[0] = i;
					ij_this_pair[1] = j;
					deletionSize_rep = tmpSize;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
	for(int i=0; i<this.pair.getDirectionalRangeList().size();i++){
	    DirectionalRange rangeI = this.pair.getDirectionalRangeList().get(i);
	    if(rangeI.isFwd()){
		for(int j=0; j<this.directionalRangeList.size();j++){
		    DirectionalRange rangeJ = this.directionalRangeList.get(j);
		    if(!rangeJ.isFwd()){
			int tmpSize = rangeJ.getMin() - rangeI.getMax() - 1;
			boolean update = true;
			if(foundDeletion_this_pair){
			    if(tmpSize <= deletionSize_this_pair)
				update = false;
			}
			if(foundDeletion_pair_this){
			    if(tmpSize <= deletionSize_pair_this)
				update = false;
			}
			if(update){
			    if( (rangeI.getMax() < rangeJ.getMin())
				&& rangeI.isSimpleDeletionDirectional(rangeJ) ){
				int lackingCoverageVal = t.isPathLackingCoverage(new IntervalPos(t, rangeI.getMax()), 
										 new IntervalPos(t, rangeJ.getMin()));
				if(lackingCoverageVal == 1){
				    foundDeletion_pair_this = true;
				    foundDeletion_this_pair = false;
				    deletionSize_pair_this = tmpSize;
				    deletionSize_this_pair = 0;
				    ij_pair_this[0] = i;
				    ij_pair_this[1] = j;
				    if(foundDeletion_this_pair_rep){
					this.checkRepetitiveDeletionAtTheEnd = false;
					this.pair.setCheckRepDeletionAgain(false);
					foundDeletion_this_pair_rep = false;
				    }
				    if(foundDeletion_pair_this_rep){
					this.checkRepetitiveDeletionAtTheEnd = false;
					this.pair.setCheckRepDeletionAgain(false);
					foundDeletion_pair_this_rep = false;
				    }
				}else if(!foundDeletion_this_pair && !foundDeletion_pair_this && lackingCoverageVal == -1){
				    if(tmpSize < deletionSize_rep){
					this.checkRepetitiveDeletionAtTheEnd = true;
					foundDeletion_pair_this_rep = true;
					foundDeletion_this_pair_rep = false;
					this.pair.setCheckRepDeletionAgain(true);
					ij_pair_this[0] = i;
					ij_pair_this[1] = j;
					deletionSize_rep = tmpSize;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
	
	int[] tmp = new int[3];
	tmp[0] = -3;
	if(foundDeletion_this_pair){
	    this.checkRepetitiveDeletionAtTheEnd = false;
	    this.pair.setCheckRepDeletionAgain(false);
	    tmp[0] = 0;
	    tmp[1] = ij_this_pair[0];
	    tmp[2] = ij_this_pair[1];
	}else if(foundDeletion_pair_this){
	    this.checkRepetitiveDeletionAtTheEnd = false;
	    this.pair.setCheckRepDeletionAgain(false);
	    tmp[0] = 1;
	    tmp[1] = ij_pair_this[0];
	    tmp[2] = ij_pair_this[1];
	}else if(this.checkRepetitiveDeletionAtTheEnd){
	    if(foundDeletion_this_pair_rep){
		tmp[0] = -1;
		tmp[1] = ij_this_pair[0];
		tmp[2] = ij_this_pair[1];
	    }else if(foundDeletion_pair_this_rep){
		tmp[0] = -2;
		tmp[1] = ij_pair_this[0];
		tmp[2] = ij_pair_this[1];
	    }
	}
	return tmp;
    }



    //return array size of 3. [0]: 0 if this-pair tandem duplication, 1 if pair-this tandem duplication, -1 if no Tandem Duplication
    public int[] isTandemDuplication(){
	boolean foundTandemDUP_this_pair = false;
	boolean foundTandemDUP_pair_this = false;

	int[] ij_this_pair = new int[2];
	int[] ij_pair_this = new int[2];
	
	int tandemDUPSize_this_pair = 0;
	int tandemDUPSize_pair_this = 0;
	
	for(int i=0; i<this.directionalRangeList.size();i++){
	    DirectionalRange rangeI = this.directionalRangeList.get(i);
	    if(!rangeI.isFwd()){ //rangeI direction has to be REV 
		for(int j=0;j<this.pair.getDirectionalRangeList().size();j++){
		    DirectionalRange rangeJ = this.pair.getDirectionalRangeList().get(j);
		    if(rangeJ.isFwd()){//rangeJ direction has to be FWD
			int tmpSize = rangeJ.getMax() - rangeI.getMin() - 1;
			boolean update = true;
			if(foundTandemDUP_this_pair){//if there was another signal, 
			    if(tmpSize >= tandemDUPSize_this_pair || tmpSize > Constants.MAX_TANDUP_SIZE)//if tmpSize is larger, we don't proceed
				update = false;
			}
			if(update){//if tmpSize is smaller, we proceed to check
			    if( (rangeJ.getMax() > rangeI.getMin())
				&& rangeI.isTandemDUPDirectional(rangeJ) ){
				foundTandemDUP_this_pair = true;
				tandemDUPSize_this_pair = tmpSize;
				ij_this_pair[0] = i;
				ij_this_pair[1] = j;
			    }
			}
		    }
		}
	    }
	}
	
	for(int i=0; i<this.pair.getDirectionalRangeList().size();i++){
	    DirectionalRange rangeI = this.pair.getDirectionalRangeList().get(i);
	    if(!rangeI.isFwd()){
		for(int j=0; j<this.directionalRangeList.size();j++){
		    DirectionalRange rangeJ = this.directionalRangeList.get(j);
		    if(rangeJ.isFwd()){
			int tmpSize = rangeJ.getMax() - rangeI.getMin() - 1;
			boolean update = true;
			if(foundTandemDUP_this_pair){
			    if(tmpSize >=tandemDUPSize_this_pair || tmpSize > Constants.MAX_TANDUP_SIZE)
				update = false;
			}
			if(foundTandemDUP_pair_this){
			    if(tmpSize >= tandemDUPSize_pair_this || tmpSize > Constants.MAX_TANDUP_SIZE)
				update = false;
			}
			if(update){
			    if( (rangeJ.getMax() > rangeI.getMin())
				&& rangeI.isTandemDUPDirectional(rangeJ) ){
				foundTandemDUP_pair_this = true;
				foundTandemDUP_this_pair = false;
				tandemDUPSize_pair_this = tmpSize;
				ij_pair_this[0] = i;
				ij_pair_this[1] = j;
			    }
			}
		    }
		}
	    }
	}
	
	int[] tmp = new int[3];
	tmp[0] = -1;
	if(foundTandemDUP_this_pair){
	    if(tandemDUPSize_this_pair >= Constants.MIN_TANDUP_SIZE){
		tmp[0] = 0;
		tmp[1] = ij_this_pair[0];
		tmp[2] = ij_this_pair[1];
	    }else
		tmp[0] = -2;
	}else if(foundTandemDUP_pair_this){
	    if(tandemDUPSize_pair_this >= Constants.MIN_TANDUP_SIZE){
		tmp[0] = 1;
		tmp[1] = ij_pair_this[0];
		tmp[2] = ij_pair_this[1];
	    }else
		tmp[0] = -2;
	}
	return tmp;
	
    }

    /* only checking for non-repetative anchoring edge deletion 
     * allow one-side to be repetative (deletion via homologous recombination)
     * this halfCluster must be on unique path.
     * returns null if NOT a simple deletion
     **/
    /*
    public ArrayList<int[]> isSimpleDeletion(Thread t){
	ArrayList<int[]> listOfIndexPairForDirectionalRanges = new ArrayList<int[]>();
	
	
	DirectionalRange anchorRange = this.directionalRangeList.get(0);

	// ----anchorRange----->gp1          gp2 <------pair HCL's DR----
	if(anchorRange.isFwd()){
	    //for each DR from pair HCL
	    //we search backaward since we are looking for largest deletion supported by coverage.
	    //Also this way we can stop as soon as two clusters meet
	    for(int j=this.pair.getDirectionalRangeList().size()-1;j>=0;j--){
		if(this.pair.getDirectionalRangeList().get(j).getMin() < anchorRange.getMax())
		    break;//we don't need to check further
		
		//direction must be opposite to anchorRange
		//and check if it has proper deletion signal
		if(!this.pair.getDirectionalRangeList().get(j).isFwd()
		   && anchorRange.isSimpleDeletion(this.pair.getDirectionalRangeList().get(j))){
		    //needs to check depth
		    IntervalPos anchor_ip_tip = new IntervalPos(t, anchorRange.getMax());
		    IntervalPos pair_ip_tip = new IntervalPos(t, this.pair.getDirectionalRangeList().get(j).getMin());
		    
		    if(t.isPathLackingCoverage(anchor_ip_tip, pair_ip_tip)){
			listOfIndexPairForDirectionalRanges.add(new int[] {0,j});
			break;
		    }
		}
	    }
	}
	// ------pair HCL's DR---->gp1          gp2<----anchorRange-----
	else{
	    //for each DR from pair HCL
	    //we search forward since we are looking for largest deletion supported by coverage.
	    //Also this way we can stop as soon as two clusters meet
	    for(int j=0;j<this.pair.getDirectionalRangeList().size();j++){
		if(this.pair.getDirectionalRangeList().get(j).getMax() > anchorRange.getMin())
		    break;//we don't need to check further
		
		//direction must be opposite to anchorRange
		//and check if it has proper deletion signal
		if(this.pair.getDirectionalRangeList().get(j).isFwd()
		   && this.pair.getDirectionalRangeList().get(j).isSimpleDeletion(anchorRange)){
		    //needs to check depth
		    IntervalPos pair_ip_tip = new IntervalPos(t, this.pair.getDirectionalRangeList().get(j).getMax());
		    IntervalPos anchor_ip_tip = new IntervalPos(t, anchorRange.getMin());
		    
		    if(t.isPathLackingCoverage(pair_ip_tip, anchor_ip_tip)){
			listOfIndexPairForDirectionalRanges.add(new int[] {0,j});
			break;
		    }
		}
	    }
	
	}
	if(listOfIndexPairForDirectionalRanges.size() > 0)
	    return listOfIndexPairForDirectionalRanges;
	return null;
    }
    */

    public String printAsDR(){
	StringBuffer bf = new StringBuffer();
	bf.append("############### HCL (" + this.id + ")####################\n");
	//System.err.println(" Multiplicity:\t" + this.directionalRangeList.size() );
	//System.err.println(" #Edges:\t" + this.numEdges());
	//System.err.println(" #Reads:\t" + this.numReads() + "\n");
	bf.append("#[E\tM\tR]:\t"+ this.numEdges() + "\t" + this.directionalRangeList.size() + "\t" + this.numReads() + "\n");
	for(int i=0;i<this.directionalRangeList.size();i++){
	    bf.append("#\tDR" + i + "\t" + this.directionalRangeList.get(i).toString()+"\n");
	}
	bf.append("#----------------HCL.pair ("+ this.pair.getId() +")--------------------\n");
	//System.err.println(" Multiplicity:\t" + this.pair.getDirectionalRangeList().size() );
	//System.err.println(" #Edges:\t" + this.pair.numEdges());
	//System.err.println(" #Reads:\t" + this.pair.numReads() + "\n");
	bf.append("#[E\tM\tR]:\t"+ this.pair.numEdges() + "\t" + this.pair.getDirectionalRangeList().size() + "\t" + this.pair.numReads() + "\n");
	for(int i=0;i<this.pair.getDirectionalRangeList().size();i++){
	    bf.append("#\tDR" + i + "\t" + this.pair.getDirectionalRangeList().get(i).toString() + "\n");
	}
	bf.append("#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n\n");
	return bf.toString();
    }


    public String print(){
	return this.print(false, null);
    }
    
    public void printDirect(){
	System.err.println("##########################################\n");
	String thisResult = this.print();
	System.err.println("------------------------------------------");
	String pairResult = this.pair.print();
	//System.out.println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
	//String pairResultSorted = this.halfClusterTable.get(curInt).getPair().print(true, this.t);//true for sorting EdgeList and t is needed to sort
	System.err.println("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
	System.err.println("U>" + thisResult + "\t" + pairResult);
    }


    public String print(boolean sortBefore, Thread t){
	if(sortBefore){
	    if(this.edgeList.size() > 1)
		this.sortEdgeList(t);
	}
	//System.out.println("ClusterID:\t" + this.id);
	StringBuffer result = new StringBuffer();
	boolean[] inclusionArr = this.getInclusionArrayForSignificantEdges();
	int numEdges = 0;
	for(int i=0;i<inclusionArr.length;i++){
	    if(inclusionArr[i])
		numEdges++;
	}
	if(numEdges == 0)
	    return null;
	System.err.println("ClusterID:\t" + this.id + "\t[numEdges : " + this.numEdges() + " ]\t[numReads : " + this.numReads() + " ]");
	result.append(this.id + "\t" + numEdges + "\t" + this.numReads(inclusionArr) );
	
	for(int i=0;i<this.edgeList.size();i++){
	    if(inclusionArr[i]){
		EdgeCluster eclu = this.edgeClusterHash.get(this.edgeList.get(i));
		//System.out.println("Edge:" + i + "\t" + this.edgeList.get(i).getId() + "\t" + eclu.minPosRead().getGraphPos() + "\t" + eclu.maxPosRead().getGraphPos());
		System.err.println("Edge[" + i + "]:\t" + this.edgeList.get(i).getId() + "\t" + eclu.minPosRead().getGraphPos() + "\t" + eclu.maxPosRead().getGraphPos() );
		System.err.println("Linear Positions :\n");
		System.err.print(this.edgeList.get(i).possibleLinearPositionPairs(eclu.minPosRead(), eclu.maxPosRead(), positionRead.isFwdAtEc(0)));
		
		result.append(this.edgeList.get(i).getFirstLinearPositionPairs(eclu.minPosRead(), eclu.maxPosRead(), positionRead.isFwdAtEc(0)));
		
	    }
	}
	System.err.println();
	
	return result.toString();
    }

    
    public String printPlain(){
	StringBuffer result = new StringBuffer();
	int numEdges = this.numEdges();
	System.err.println("ClusterID:\t" + this.id + "\t[numEdges : " + numEdges + " ]\t[numReads : " + this.numReads() + " ]");
	result.append(this.id + "\t" + numEdges + "\t" + this.numReads() );
	for(int i=0;i<this.edgeList.size();i++){
	    //if(inclusionArr[i]){
	    EdgeCluster eclu = this.edgeClusterHash.get(this.edgeList.get(i));
	    //System.out.println("Edge:" + i + "\t" + this.edgeList.get(i).getId() + "\t" + eclu.minPosRead().getGraphPos() + "\t" + eclu.maxPosRead().getGraphPos());
	    System.err.println("Edge[" + i + "]:\t" + this.edgeList.get(i).getId() + "\t" + eclu.minPosRead().getGraphPos() + "\t" + eclu.maxPosRead().getGraphPos() );
	    System.err.println("Linear Positions :\n");
	    System.err.print(this.edgeList.get(i).possibleLinearPositionPairs(eclu.minPosRead(), eclu.maxPosRead(), positionRead.isFwdAtEc(0)));
	    
	    result.append(this.edgeList.get(i).getFirstLinearPositionPairs(eclu.minPosRead(), eclu.maxPosRead(), positionRead.isFwdAtEc(0)));
	    
		//}
	}
	System.err.println();
	
	return result.toString();
    
    }



    /* public void print(){
	System.out.println(this.getVitalString() + "\t" + this.getPair().getVitalString() + "\t")
    }

    public String getVitalString(){
	return this.id + "\t#Edges:"+ this.numEdges() + "\t#Reads:" + this.numReads();
	}*/

    //STILL UPDATING: WORK IN PROGRESS. FOR PRINTING CLUSTER
    public String[] getTokenizedString(){
	//0: clusterID
	//1: numEdges
	//2: numReads
	//3: edges
	String[] results = new String[10];
	return results;
    }

    public int numEdges(){
	return this.edgeList.size();
    }
    
    public int numReads(){
	int size = 0;
	
	for(int i=0;i<edgeList.size();i++)
	    size += this.edgeClusterHash.get(edgeList.get(i)).numReads();
	
	return size;
    }

    public int numReads(boolean[] inclusionArr){
	int size = 0;
	for(int i=0;i<edgeList.size();i++)
	    if(inclusionArr[i])
		size += this.edgeClusterHash.get(edgeList.get(i)).numReads();
	
	return size;
    }
    
    private static int nextId = 1;


    private boolean checkRepetitiveDeletionAtTheEnd;
    
    private int id;

    private int curLinPos;
    
    private HalfCluster pair; // this is the pairing half-cluster
    private Read positionRead; 

    private ArrayList<Edge> edgeList; // this contains all the Edges that are spanned by this halfcluster
    private HashMap<Edge, EdgeCluster> edgeClusterHash; //Fetech EdgeCluster for each Edge.

    private boolean leading;
    
    private ArrayList<DirectionalRange> directionalRangeList;

    private int tipIndex; //index for accessing tip Edge in edgeList.
    private int stemIndex; //index for accessing stem Edge in edgeList.


        public String getEventTypeOLD(HalfCluster other, boolean debug){
	int tDRsize = this.directionalRangeList.size();
	int tpDRsize = this.pair.getDirectionalRangeList().size();
	int oDRsize = other.getDirectionalRangeList().size();
	int opDRsize = other.getPair().getDirectionalRangeList().size();
	
	if(debug){
	    System.err.println("=T======================================");
	    this.printPlain();
	    System.err.println("=TP======================================");	
	    this.pair.printPlain();
	    System.err.println("=O======================================");
	    other.printPlain();
	    System.err.println("=OP======================================");
	    other.getPair().printPlain();
	    System.err.println("=4======================================");
	}
	if(debug){
	    System.err.println("tDRsize = " + tDRsize );
	    System.err.println("tpDRsize = " + tDRsize);
	    System.err.println("oDRsize = " + oDRsize);
	    System.err.println("opDRsize = " + opDRsize + "\n");
	}
	 

	int[] tRelativeRelationVal = new int[tDRsize*tpDRsize];
	int[] oRelativeRelationVal = new int[oDRsize*opDRsize];
	
	for(int i=0; i<tDRsize; i++){
	    for(int j=0; j<tpDRsize; j++){
		tRelativeRelationVal[i*tpDRsize + j] = this.getNthDirectionalRange(i).getRelativeRelationType(this.pair.getNthDirectionalRange(j));
		if(debug){
		    System.err.println("getting relativeRelation for\ntDR :\t" + this.getNthDirectionalRange(i).toString() + "\ntpDR :\t" + this.pair.getNthDirectionalRange(j).toString());
		    System.err.println("RelationTYPE:\t" + tRelativeRelationVal[i*tpDRsize + j]);
		    System.err.println();
		}
	    }
	}
	
	

	for(int i=0; i<oDRsize; i++){
	    for(int j=0; j<opDRsize; j++){
		if(debug)
		  System.err.println("i:j\t" + i + "\t" + j);
		oRelativeRelationVal[i*opDRsize + j] = other.getNthDirectionalRange(i).getRelativeRelationType(other.getPair().getNthDirectionalRange(j));
		if(debug){
		    System.err.println("getting relativeRelation for\noDR :\t" + other.getNthDirectionalRange(i).toString() + "\nopDR :\t" + other.getPair().getNthDirectionalRange(j).toString());
		    System.err.println("RelationTYPE:\t" + oRelativeRelationVal[i*opDRsize + j]);
		}
	    }
	}

	
	//from i, accessing this DR: ab
	for(int i=0; i<tRelativeRelationVal.length; i++){
	    DirectionalRange tDR = this.getNthDirectionalRange(i/tpDRsize);
	    DirectionalRange tpDR = this.getPair().getNthDirectionalRange(i%tpDRsize);
	    for(int j=0; j<oRelativeRelationVal.length; j++){
		//DirectionalRange = tDR = this.getNthDirectionalRange(i/tpDRsize);
		//DirectionalRange = tpDR = this.getPair.getNthDirectionalRange(i%tpDRsize);
		DirectionalRange oDR = other.getNthDirectionalRange(j/opDRsize);
		DirectionalRange opDR = other.getPair().getNthDirectionalRange(j%opDRsize);

		
		//#3 Directed Transposition
		if( (tRelativeRelationVal[i] == 4       //[NO, TP, +-] --this--> <--pair--
		     || tRelativeRelationVal[i] == 12   //[OV, TP, +-] -<P=this,pair=T>-
		     || tRelativeRelationVal[i] == 15   //[OV, PT, -+] -<P=pair,this=T>-
		     || tRelativeRelationVal[i] == 7    //[NO, PT, -+] <--pair-- --this--> 
		     ) &&
		    (oRelativeRelationVal[j] == 5      //[NO, OP, -+] <--other-- --pair-->
		     || oRelativeRelationVal[j] == 13  //[OV, OP, -+] -<O=other,pair=P>-
		     || oRelativeRelationVal[j] == 14  //[OV, PO, +-] -<O=pair,other=P>-
		     || oRelativeRelationVal[j] == 6   //[NO, PO, +-] --pair--> <--other--
		     ))		    
		    {
			//
			// <--other---   ----this-->
			// -------IS-ELEMENT--------
			//
			// directed insertion 
			//
			// ---other.p--> <---this.p---
			// -------------|-------------
			// insertion site "|"
			//
			if(oDR.isBoundaryDefinerFwd(tDR) && opDR.isWithinFacingFwd(tpDR)){
			    System.out.println("TRANSPOSTION\t"+"3-1"+"\t" + other.getNumReadsDRSize() + this.getNumReadsDRSize() 
					       + oDR.toString() +"\t"+ tDR.toString() + "\t|\t" + opDR.toString()+"\t" + tpDR.toString()
					       + "\t" + other.getPair().getNumReadsDRSize() + this.pair.getNumReadsDRSize());
			    return "3-1";//ASSIGN_TRANSPOSITION_3-1;
			}
			//
			// <--this.p---   ----other.p-->
			// -------IS-ELEMENT--------
			//
			// directed insertion 
			//
			// ---this---> <---other---
			// -----------|-------------
			// insertion site "\t|\t"
			//
			else if(tpDR.isBoundaryDefinerFwd(opDR) && tDR.isWithinFacingFwd(oDR)){
			    System.out.println("TRANSPOSTION\t"+"3-2"+"\t" + this.pair.getNumReadsDRSize() + other.getPair().getNumReadsDRSize() 
					       + tpDR.toString() +"\t"+ opDR.toString() + "\t|\t" + tDR.toString() +"\t"+ oDR.toString() 
					       + "\t" + this.getNumReadsDRSize() +other.getNumReadsDRSize());
			    return "3-2";//ASSIGN_TRANSPOSITION_3-2;
			}
			
		    }
		//#1 Directed Transposition
		else if( (tRelativeRelationVal[i] == 5      //[NO, TP, -+] <--this-- --pair-->
			  || tRelativeRelationVal[i] == 13  //[OV, TP, -+] -<T=this,pair=P>-
			  || tRelativeRelationVal[i] == 14  //[OV, PT, +-] -<T=pair,this=P>-
			  || tRelativeRelationVal[i] == 6   //[NO, PT, +-] --pair--> <--this--
			  ) &&
			 (oRelativeRelationVal[j] == 4       //[NO, TP, +-] --this--> <--pair--
			  || oRelativeRelationVal[j] == 12   //[OV, TP, +-] -<P=this,pair=T>-
			  || oRelativeRelationVal[j] == 15   //[OV, PT, -+] -<P=pair,this=T>-
			  || oRelativeRelationVal[j] == 7    //[NO, PT, -+] <--pair-- --this--> 
			  ))
		    {
			//
			// <--this---   ----other-->
			// -------IS-ELEMENT--------
			//
			// directed insertion 
			//
			// ---this.p--> <---other.p---
			// ------------|-------------
			// insertion site "\t|\t"
			//
			if(tDR.isBoundaryDefinerFwd(oDR) && tpDR.isWithinFacingFwd(opDR)){
			    System.out.println("TRANSPOSTION\t"+"1-1"+"\t" + this.getNumReadsDRSize() +other.getNumReadsDRSize()
					       + tDR.toString() +"\t"+ oDR.toString() + "\t|\t" + tpDR.toString() +"\t"+ opDR.toString()
					       + "\t" + this.pair.getNumReadsDRSize() + other.getPair().getNumReadsDRSize());
			    return "1-1";//ASSIGN_TRANSPOSITION_1-1;
			}
			//
			// <--other.p---   ----this.p-->
			// -------IS-ELEMENT--------
			//
			// directed insertion 
			//
			// ---other---> <---this---
			// ------------|-------------
			// insertion site "\t|\t"
			//
			else if(opDR.isBoundaryDefinerFwd(tpDR) && oDR.isWithinFacingFwd(tDR)){
			    System.out.println("TRANSPOSTION\t"+"1-2"+"\t" + other.getPair().getNumReadsDRSize() + this.pair.getNumReadsDRSize()
					       + opDR.toString() +"\t"+ tpDR.toString() + "\t|\t" + oDR.toString() +"\t"+ tDR.toString()
					       + "\t" + other.getNumReadsDRSize() + this.getNumReadsDRSize());
			    return "1-2";//ASSIGN_TRANSPOSITION_1-2;
			}
		    }
		//#2 Directed Transposition
		else if( (tRelativeRelationVal[i] == 5      //[NO, TP, -+] <--this-- --pair-->
			  || tRelativeRelationVal[i] == 13  //[OV, TP, -+] -<T=this,pair=P>-
			  || tRelativeRelationVal[i] == 14  //[OV, PT, +-] -<T=pair,this=P>-
			  || tRelativeRelationVal[i] == 6   //[NO, PT, +-] --pair--> <--this--
			  ) &&
			 (oRelativeRelationVal[j] == 5       //[NO, TP, +-] --this--> <--pair--
			  || oRelativeRelationVal[j] == 13   //[OV, TP, +-] -<P=this,pair=T>-
			  || oRelativeRelationVal[j] == 14   //[OV, PT, -+] -<P=pair,this=T>-
			  || oRelativeRelationVal[j] == 6    //[NO, PT, -+] <--pair-- --this--> 
			  ))
		    {
			//
			// <--this---   ----other.p-->
			// -------IS-ELEMENT--------
			//
			// directed insertion 
			//
			// ---this.p--> <---other----
			// ------------|-------------
			// insertion site "\t|\t"
			//
			if(tDR.isBoundaryDefinerFwd(opDR) && tpDR.isWithinFacingFwd(oDR)){
			    System.out.println("TRANSPOSTION\t"+"2-1"+"\t" + this.getNumReadsDRSize() + other.getPair().getNumReadsDRSize()
					       + tDR.toString() +"\t"+ opDR.toString() + "\t|\t" + tpDR.toString() +"\t"+ oDR.toString()
					       + "\t" + this.pair.getNumReadsDRSize() + other.getNumReadsDRSize());
			    return "2-1";//ASSIGN_TRANSPOSITION_2-1;
			}
			//
			// <--other---   ----this.p-->
			// -------IS-ELEMENT--------
			//
			// directed insertion 
			//
			// ---other.p--> <---this----
			// ------------|-------------
			// insertion site "\t|\t"
			//
			else if(oDR.isBoundaryDefinerFwd(tpDR) && opDR.isWithinFacingFwd(tDR)){
			    System.out.println("TRANSPOSTION\t"+"2-2"+"\t" + other.getNumReadsDRSize() + this.pair.getNumReadsDRSize() 
					       + oDR.toString() +"\t"+ tpDR.toString() + "\t|\t" + opDR.toString() +"\t"+ tDR.toString()
					       + "\t" + other.getPair().getNumReadsDRSize() + this.getNumReadsDRSize());
			    return "2-2";//ASSIGN_TRANSPOSITION_2-2;
			}
		    }
		//#4 Directed Transposition
		else if( (tRelativeRelationVal[i] == 4      //[NO, TP, -+] <--this-- --pair-->
			  || tRelativeRelationVal[i] == 12  //[OV, TP, -+] -<T=this,pair=P>-
			  || tRelativeRelationVal[i] == 15  //[OV, PT, +-] -<T=pair,this=P>-
			  || tRelativeRelationVal[i] == 7   //[NO, PT, +-] --pair--> <--this--
			  ) &&
			 (oRelativeRelationVal[j] == 4       //[NO, TP, -+] <--this-- --pair-->
			  || oRelativeRelationVal[j] == 12  //[OV, TP, -+] -<T=this,pair=P>-
			  || oRelativeRelationVal[j] == 15  //[OV, PT, +-] -<T=pair,this=P>-
			  || oRelativeRelationVal[j] == 7   //[NO, PT, +-] --pair--> <--this--
			  ))
		    {
			//
			// <--this.p---   ----other-->
			// -------IS-ELEMENT--------
			//
			// directed insertion 
			//
			// ---this----> <---other.p---
			// ------------|-------------
			// insertion site "\t|\t"
			//
			if(tpDR.isBoundaryDefinerFwd(oDR) && tDR.isWithinFacingFwd(oDR)){
			    System.out.println("TRANSPOSTION\t"+"4-1"+"\t" + this.pair.getNumReadsDRSize() + other.getNumReadsDRSize()
					       + tpDR.toString() +"\t"+ oDR.toString() + "\t|\t" + tDR.toString() +"\t"+ opDR.toString()
					       + "\t" + this.getNumReadsDRSize() + other.getPair().getNumReadsDRSize());
			    return "4-1";//ASSIGN_TRANSPOSITION_4-1;
			}
			//
			// <--other.p---   ----this-->
			// -------IS-ELEMENT--------
			//
			// directed insertion 
			//
			// ---other----> <---this.p---
			// -------------|-------------
			// insertion site "\t|\t"
			//
			else if(opDR.isBoundaryDefinerFwd(tDR) && oDR.isWithinFacingFwd(tpDR)){
			    System.out.println("TRANSPOSTION\t"+"4-2"+"\t" + other.getPair().getNumReadsDRSize() + this.getNumReadsDRSize()
					       + opDR.toString() +"\t"+ tDR.toString() + "\t|\t" + oDR.toString() +"\t"+ tpDR.toString()
					       + "\t" + other.getNumReadsDRSize() + this.pair.getNumReadsDRSize());
			    return "4-2";//ASSIGN_TRANSPOSITION_4-2;
			}
		    }
		//#5 inverted transposition && INVERSION
		else if( (tRelativeRelationVal[i] == 1      //[NO, TP, --] <--this-- <--pair--
			  || tRelativeRelationVal[i] == 9   //[OV, TP, --] <<=this,pair=-
			  || tRelativeRelationVal[i] == 10  //[OV, PT, --] <<=pair,this=- 
			  || tRelativeRelationVal[i] == 2   //[NO, PT, --] <--pair-- <--this--
			  ) &&
			 (oRelativeRelationVal[j] == 0       //[NO, TP, ++] --this--> --pair-->
			  || oRelativeRelationVal[j] == 8    //[OV, TP, ++] -=this,pair=>>
			  || oRelativeRelationVal[j] == 11   //[OV, PT, ++] -=pair,this=>>
			  || oRelativeRelationVal[j] == 3    //[NO, PT, ++] --pair--> --this-->
			  ))
		    {
		
			//
			// <--this---   ----other-->
			// -------IS-ELEMENT--------
			//
			// directed insertion 
			//
			// ---other.p--> <---this.p---
			// -------------|-------------
			// insertion site "\t|\t"
			//
			if(tDR.isBoundaryDefinerFwd(oDR) && opDR.isWithinFacingFwd(tpDR)){
			    System.out.println("TRANSPOSTION(I)\t"+"5-1"+"\t" + this.getNumReadsDRSize() + other.getNumReadsDRSize()
					       + tDR.toString() +"\t" + oDR.toString() + "\t|\t" + opDR.toString()+"\t" + tpDR.toString()
					       + "\t" + other.getPair().getNumReadsDRSize() + this.pair.getNumReadsDRSize());
			    return "5-1";//ASSIGN_TRANSPOSITION_INVERSION_5-1;
			}
			//
			// <--this---   ----other.p-->
			// -------IS-ELEMENT--------
			//
			// directed insertion 
			//
			// ---other--> <---this.p---
			// -----------|-------------
			// insertion site "\t|\t"
			//
			else if(tDR.isBoundaryDefinerFwd(opDR) && oDR.isWithinFacingFwd(tpDR)){
			    System.out.println("TRANSPOSTION(I)\t"+"5-2"+"\t" + this.getNumReadsDRSize() + other.getPair().getNumReadsDRSize() 
					       + tDR.toString() + "\t" +opDR.toString() + "\t|\t" + oDR.toString() +"\t"+ tpDR.toString()
					       + "\t" + other.getNumReadsDRSize() + this.pair.getNumReadsDRSize());
			    return "5-2";//ASSIGN_TRANSPOSITION_INVERSION_5-2;
			}
			//
			// <--this.p---   ----other-->
			// -------IS-ELEMENT--------
			//
			// directed insertion 
			//
			// ---other.p--> <---this---
			// -------------|-------------
			// insertion site "\t|\t"
			//
			else if(tpDR.isBoundaryDefinerFwd(oDR) && opDR.isWithinFacingFwd(tDR)){
			    System.out.println("TRANSPOSTION(I)\t"+"5-3"+"\t" + this.pair.getNumReadsDRSize() + other.getNumReadsDRSize()
					       + tpDR.toString() + "\t" +oDR.toString() + "\t|\t" + opDR.toString() +"\t"+ tDR.toString()
					       + "\t" + other.getPair().getNumReadsDRSize() + this.getNumReadsDRSize());
			    return "5-3";//ASSIGN_TRANSPOSITION_INVERSION_5-3;
			}
			//
			// <--this.p---   ----other.p-->
			// -------IS-ELEMENT--------
			//
			// directed insertion 
			//
			// ---other---> <---this---
			// ------------|-------------
			// insertion site "\t|\t"
			//
			else if(tpDR.isBoundaryDefinerFwd(opDR) && oDR.isWithinFacingFwd(tDR)){
			    System.out.println("TRANSPOSTION(I)\t"+"5-4"+"\t" + this.pair.getNumReadsDRSize() + other.getPair().getNumReadsDRSize()
					       + tpDR.toString() +"\t"+ opDR.toString() + "\t|\t" + oDR.toString() +"\t"+ tDR.toString()
					       + "\t" + other.getNumReadsDRSize() + this.getNumReadsDRSize() );
			    return "5-4";//ASSIGN_TRANSPOSITION_INVERSION_5-4;
			}
			//
			// INVERSIONS
			//
			//
			// --other-->     ---other.p--->
			//           <---this-----       <--this.p--
			//          |-----INVERSION-----|
			else if( (tRelativeRelationVal[i] == 1
				  || tRelativeRelationVal[i] == 9 
				   ) &&
				 (oRelativeRelationVal[j] == 0
				  || oRelativeRelationVal[j] == 8)
				 ){
			    //if(debug){
			    //		System.err.println("\t\t#### IN INVERSION 5-A");
			    //}
			    int padding = oDR.getPaddingWhenTestingisWithinFacingFwd(tDR, opDR, tpDR, true);
			    //if(debug){
			    //		System.err.println("\t\t#### IN INVERSION 5-A");
			    //		System.err.println("\\t\t\tPADDING SIZE");
			    //}
			    if(oDR.isWithinFacingFwdINVERSION(tDR, padding) && opDR.isWithinFacingFwdINVERSION(tpDR, padding)){
				System.out.println("INVERSION\t"+"5-A"+"\t" + other.getNumReadsDRSize() + this.getNumReadsDRSize()
						   + oDR.toString() +"\t"+ tDR.toString() + "\t|\t" + opDR.toString() +"\t"+ tpDR.toString()
						   + "\t" + other.getPair().getNumReadsDRSize() + this.pair.getNumReadsDRSize());
				return "5-A";//ASSIGN_INVERSION_5_A;
			    }
			}
			// --other-->     ---other.p--->
			//           <---this.p----      <--this---
			//          |-----INVERSION-----|
			else if( (tRelativeRelationVal[i] == 2
				  || tRelativeRelationVal[i] == 10 
				  ) &&
				 (oRelativeRelationVal[j] == 0
				  || oRelativeRelationVal[j] == 8)
				 ){
			    int padding = oDR.getPaddingWhenTestingisWithinFacingFwd(tpDR, opDR, tDR, true);
			    if(oDR.isWithinFacingFwdINVERSION(tpDR, padding) && opDR.isWithinFacingFwdINVERSION(tDR, padding)){
				System.out.println("INVERSION\t"+"5-B"+"\t" + other.getNumReadsDRSize() + this.pair.getNumReadsDRSize()
						   + oDR.toString() +"\t"+ tpDR.toString() + "\t|\t" + opDR.toString() +"\t"+ tDR.toString()
						   + "\t" + other.getPair().getNumReadsDRSize() + this.getNumReadsDRSize());
				return "5-B";//ASSIGN_INVERSION_5_B;
			    }
			}
			// --other.p-->     ---other----->
			//             <----this-----      <--this.p---
			//            |-----INVERSION-----|
			else if( (tRelativeRelationVal[i] == 1
				  || tRelativeRelationVal[i] == 9
				   ) &&
				 (oRelativeRelationVal[j] == 3
				  || oRelativeRelationVal[j] == 11)
				  ){
			    int padding = opDR.getPaddingWhenTestingisWithinFacingFwd(tDR, oDR, tpDR, true);
			    if(opDR.isWithinFacingFwdINVERSION(tDR, padding) && oDR.isWithinFacingFwdINVERSION(tpDR, padding)){
				System.out.println("INVERSION\t"+"5-C"+"\t" + other.getPair().getNumReadsDRSize() + this.getNumReadsDRSize()
						   + opDR.toString() +"\t"+ tDR.toString() + "\t|\t" + oDR.toString() +"\t"+ tpDR.toString()
						   + "\t" + other.getNumReadsDRSize() + this.pair.getNumReadsDRSize() );
				return "5-C";//ASSIGN_INVERSION_5_C;
			    }
			}
			// --other.p-->     ---other----->
			//             <----this.p---      <--this---
			//            |-----INVERSION-----|
			else if( (tRelativeRelationVal[i] == 2
				  || tRelativeRelationVal[i] == 10
				   ) &&
				 (oRelativeRelationVal[j] == 3
				  || oRelativeRelationVal[j] == 11)
				  ){
			    int padding = opDR.getPaddingWhenTestingisWithinFacingFwd(tpDR, oDR, tDR, true);
			    if(opDR.isWithinFacingFwdINVERSION(tpDR, padding) && oDR.isWithinFacingFwdINVERSION(tDR, padding)){
				System.out.println("INVERSION\t"+"5-D"+"\t" + other.getPair().getNumReadsDRSize() + this.pair.getNumReadsDRSize()
						   + opDR.toString() +"\t"+ tpDR.toString() + "\t|\t" + oDR.toString() +"\t"+ tDR.toString()
						   + "\t" + other.getNumReadsDRSize() + this.getNumReadsDRSize());
				return "5-D";//ASSIGN_INVERSION_5_D;
			    }
			}
		    }
		
		//#6 inverted transposition && INVERSION
		else if( (tRelativeRelationVal[i] == 0       //[NO, TP, ++] --this--> --pair-->
			  || tRelativeRelationVal[i] == 8    //[OV, TP, ++] -=this,pair=>>
			  || tRelativeRelationVal[i] == 11   //[OV, PT, ++] -=pair,this=>>
			  || tRelativeRelationVal[i] == 3    //[NO, PT, ++] --pair--> --this-->
			  ) &&
			 (oRelativeRelationVal[j] == 1      //[NO, OP, --] <--other-- <--pair--
			  || oRelativeRelationVal[j] == 9   //[OV, OP, --] <<=other,pair=-
			  || oRelativeRelationVal[j] == 10  //[OV, PO, --] <<=pair,other=- 
			  || oRelativeRelationVal[j] == 2   //[NO, PO, --] <--pair-- <--other--
			  ))
		    {
			//
			// <--other---   ----this-->
			// -------IS-ELEMENT--------
			//
			// directed insertion 
			//
			// ---this.p--> <---other.p---
			// ------------|-------------
			// insertion site "\t|\t"
			//
			if(oDR.isBoundaryDefinerFwd(tDR) && tpDR.isWithinFacingFwd(opDR)){
			    System.out.println("TRANSPOSTION(I)\t"+"6-1"+"\t" + other.getNumReadsDRSize() + this.getNumReadsDRSize()
					       + oDR.toString() +"\t"+ tDR.toString() + "\t|\t" + tpDR.toString()+"\t" + opDR.toString()
					       + "\t" + this.pair.getNumReadsDRSize() + other.getPair().getNumReadsDRSize());
			    return "6-1";//ASSIGN_TRANSPOSITION_INVERTED_6-1;
			}
			//
			// <--other---   ----this.p-->
			// -------IS-ELEMENT--------
			//
			// directed insertion 
			//
			// ---this--> <---other.p---
			// ----------|-------------
			// insertion site "\t|\t"
			//
			else if(oDR.isBoundaryDefinerFwd(tpDR) && tDR.isWithinFacingFwd(opDR)){
			    System.out.println("TRANSPOSTION(I)\t"+"6-2"+"\t" +other.getNumReadsDRSize() +  this.pair.getNumReadsDRSize()
					       + oDR.toString()+"\t" + tpDR.toString() + "\t|\t" + tDR.toString() +"\t"+ opDR.toString()
					       + "\t" + this.getNumReadsDRSize()  + other.getPair().getNumReadsDRSize());
			    return "6-2";//ASSIGN_TRANSPOSITION_INVERTED_6-2;
			}
			//
			// <--other.p---   ----this-->
			// -------IS-ELEMENT--------
			//
			// directed insertion 
			//
			// ---this.p--> <---other---
			// ------------|-------------
			// insertion site "\t|\t"
			//
			else if(opDR.isBoundaryDefinerFwd(tDR) && tpDR.isWithinFacingFwd(oDR)){
			    System.out.println("TRANSPOSTION(I)\t"+"6-3"+"\t" +other.getPair().getNumReadsDRSize() + this.getNumReadsDRSize() 
					       + opDR.toString() +"\t"+ tDR.toString() + "\t|\t" + tpDR.toString() +"\t"+ oDR.toString()
					       + "\t" + this.pair.getNumReadsDRSize() + other.getNumReadsDRSize());
			    return "6-3";//ASSIGN_TRANSPOSITION_INVERTED_6-3;
			}
			//
			// <--other.p---   ----this.pair-->
			// -------IS-ELEMENT--------
			//
			// directed insertion 
			//
			// ---this--> <---other---
			// ----------|-------------
			// insertion site "\t|\t"
			//
			else if(opDR.isBoundaryDefinerFwd(tpDR) && tDR.isWithinFacingFwd(oDR)){
			    System.out.println("TRANSPOSTION(I)\t"+"6-4"+"\t" + other.getPair().getNumReadsDRSize() + this.pair.getNumReadsDRSize()
					       + opDR.toString() +"\t"+ tpDR.toString() + "\t|\t" + tDR.toString() +"\t"+ oDR.toString()
					       + "\t" + this.getNumReadsDRSize() + other.getNumReadsDRSize());
			    return "6-4";//ASSIGN_TRANSPOSITION_INVERTED_6-4;
			}
			//
			// --this-->     -----this.p--->
			//           <---other-----      <--other.p--
			//          |-----INVERSION-----|
			else if( (tRelativeRelationVal[i] == 0
				  || tRelativeRelationVal[i] == 8 
				   ) &&
				 (oRelativeRelationVal[j] == 1
				  || oRelativeRelationVal[j] == 9)
				 ){
			    int padding = tDR.getPaddingWhenTestingisWithinFacingFwd(oDR, tpDR, opDR, true);
			    //int padding2 = tpDR.getPaddingWhenTestingisWithinFacingFwd(opDR, tDR, oDR, false);
			    if(tDR.isWithinFacingFwdINVERSION(oDR, padding) && tpDR.isWithinFacingFwdINVERSION(opDR, padding)){
				System.out.println("INVERSION\t"+"6-A"+"\t" + this.getNumReadsDRSize() + other.getNumReadsDRSize()
						   + tDR.toString() +"\t"+ oDR.toString() + "\t|\t" + tpDR.toString() +"\t"+ opDR.toString()
						   + "\t" + this.pair.getNumReadsDRSize() + other.getPair().getNumReadsDRSize());
				return "6-A";//ASSIGN_INVERSION_6_A;
			    }
			}
			//
			// --this-->     -----this.p--->
			//           <---other.p----     <--other----
			//          |-----INVERSION-----|
			else if( (tRelativeRelationVal[i] == 0
				  || tRelativeRelationVal[i] == 8 
				  ) &&
				 (oRelativeRelationVal[j] == 2
				  || oRelativeRelationVal[j] == 10)
				  ){
			    int padding = tDR.getPaddingWhenTestingisWithinFacingFwd(opDR, tpDR, oDR, true);
			    //int padding2 = tpDR.getPaddingWhenTestingisWithinFacingFwd(oDR, tDR, opDR, false);
			    if(tDR.isWithinFacingFwdINVERSION(opDR,padding) && tpDR.isWithinFacingFwdINVERSION(oDR,padding)){
				System.out.println("INVERSION\t"+"6-B"+"\t" + this.getNumReadsDRSize() + other.getPair().getNumReadsDRSize()
						   + tDR.toString() +"\t"+ opDR.toString() + "\t|\t" + tpDR.toString()+"\t" + oDR.toString()
						   + "\t" + this.pair.getNumReadsDRSize() + other.getNumReadsDRSize());
				return "6-B";//ASSIGN_INVERSION_6_B;
			    }
			}
			//
			// --this.p->     -----this---->
			//           <----other----      <--other.p----
			//          |-----INVERSION-----|
			else if( (tRelativeRelationVal[i] == 3
				  || tRelativeRelationVal[i] == 11
				   ) &&
				 (oRelativeRelationVal[j] == 1
				  || oRelativeRelationVal[j] == 9)
				 ){
			    int padding = tpDR.getPaddingWhenTestingisWithinFacingFwd(oDR, tDR, opDR, true);
			    if(tpDR.isWithinFacingFwdINVERSION(oDR, padding) && tDR.isWithinFacingFwdINVERSION(opDR, padding)){
				System.out.println("INVERSION\t"+"6-C"+"\t" + this.pair.getNumReadsDRSize() + other.getNumReadsDRSize()
						   + tpDR.toString()+"\t" + oDR.toString() + "\t|\t" + tDR.toString() +"\t"+ opDR.toString()
						   + "\t" + this.getNumReadsDRSize() + other.getPair().getNumReadsDRSize());
				return "6-C";//ASSIGN_INVERSION_6_C;
			    }
			}
			//
			// --this.p->     -----this---->
			//           <----other.p----    <--other-----
			//          |-----INVERSION-----|
			else if( (tRelativeRelationVal[i] == 3
				  || tRelativeRelationVal[i] == 11 
				   ) &&
				 (oRelativeRelationVal[j] == 2
				  || oRelativeRelationVal[j] == 10)
				  ){
			    int padding = tpDR.getPaddingWhenTestingisWithinFacingFwd(opDR, tDR, oDR, true);
			    if(tpDR.isWithinFacingFwdINVERSION(opDR, padding) && tDR.isWithinFacingFwdINVERSION(oDR, padding)){
				System.out.println("INVERSION\t"+"6-D"+"\t" + this.pair.getNumReadsDRSize() + other.getPair().getNumReadsDRSize()
						   + tpDR.toString() +"\t"+ opDR.toString() + "\t|\t" + tDR.toString() +"\t"+ oDR.toString()
						   + "\t" + this.getNumReadsDRSize() + other.getNumReadsDRSize());
				return "6-D";//ASSIGN_INVERSION_6_D;
			    }
			}
			
		    }
		
		
		//if(tRelativeRelationVal[i] == 0){
		//   if(oRelativeRelationVal[j] == 1){
			//checkinf for case   ----this----> <----other---              ----this.pair ---> <---other.pair----
		//	if(this.getNthDirectionalRange(i/tpDRsize).isWithinFacingFwd(other.getNthDirectionalRagen(j/opDRsize))
		//	   && this.getNthDirectionalRange(i%tpDRsize).isWithinFacingFwd(other.getNthDirectionalRange(j%opDRsize)))
		//		    {
		//		return ASSIGNINVERSION1;
		//		    }
		//   }else if(oRelatiiveRelationVal[j] == 2){
		//	//checking for case   ----this----> <---other.pair---          ----this.pair-----> <---other----
		//	if(this.getNthDirectionalRange(i/tpDRsize).isWithinFacingFwd(other.getPair().getNthDirectionalRagen(j%opDRsize))
		//	   && this.getPair().getNthDirectionalRange(i%tpDRsize).isWithinFacingFwd(other.getNthDirectionalRange(j/opDRsize)))
		//	    {
		//		return ASSIGNINVERSION2;
		//	    }
		//  }
		//}else if(tRelativeRelationVal[j] == 3){
		//    if(oRelativeRelationVal[j] == 1){
		//		//checking for case  ----this.pair---> <---other----           ----this----> <---other.pair---
		//	if(this.getPair().getNthDirectionalRange(i%tpDRsize).isWithinFacingFwd(other.getNthDirectionalRagen(j/opDRsize))
		//	   && this.getNthDirectionalRange(i/tpDRsize).isWithinFacingFwd(other.getPair().getNthDirectionalRange(j%opDRsize)))
		//	    {
		//		return ASSIGNINVERSION3;
		//	    }
		//  }else if(oRelativeRelationVal[j] == 2){
		//	//checkinf for case ----this.pair---> <----other.pair---          ---this----> <---other---
		//	if(this.getNthDirectionalRange(i%tpDRsize).isWithinFacingFwd(other.getNthDirectionalRagen(j/opDRsize))
		//	   && this.getNthDirectionalRange(i/tpDRsize).isWithinFacingFwd(other.getNthDirectionalRange(j%opDRsize)))
		//	    {
		//		return ASSIGNINVERSION4;
		//	    }
		//  }
		//    }
	    }
	}
	return null;
    }



}

class EdgeCluster{
    
    
    /* 
     * minPosRead and maxPosRead are defined based on graphPos 
     * --> meaning if an interval direction of this edge is reverse, min is max, max is min.
     */
    public EdgeCluster(){
	this.minPosIndex = 0;
	this.maxPosIndex = 0;
	this.readList = new ArrayList<Read>();
	this.totalBases = 0;
    }


    public int numReads(){
	return this.readList.size();
    }

    public int getTotalBases(){
	return this.totalBases;
    }

    private void updateTotalBases(Read r){
	this.totalBases += r.getLength();
    }
    
    

    public void addReadsFromEdgeCluster(EdgeCluster eClu){
	//update minPosRead
	if(this.minPosRead().compareToReadOnSameEdge(eClu.minPosRead()) > 0)
	    this.minPosIndex = this.readList.size() + eClu.getMinPosIndex();
	//update maxPosRead
	if(this.maxPosRead().compareToReadOnSameEdge(eClu.maxPosRead()) < 0)
	    this.maxPosIndex = this.readList.size() + eClu.getMaxPosIndex();
	
	//now merge the list
	for(int i=0;i<eClu.getReadList().size();i++){
	    //this.readList.add(eClu.getReadList().get(i));
	    this.addReadSimple(eClu.getReadList().get(i));
	}

    }
    
    //should always use this method to addRead to readList.
    //NEVER CALL readList.add(r);
    private void addReadSimple(Read r){
	this.readList.add(r);
	this.updateTotalBases(r);
    }

    // Adding a read to Edge Cluster
    // update min and max if needed.
    public void addRead(Read r){
	//this.readList.add(r);
	this.addReadSimple(r);
	if(this.readList.size() > 1){
	    int curIndex = this.readList.size() - 1;
	    
	    //if smaller than current min, update minPos
	    if(this.readList.get(this.minPosIndex).compareToReadOnSameEdge(r) > 0)
		this.minPosIndex = curIndex;
	    else if(this.readList.get(this.maxPosIndex).compareToReadOnSameEdge(r) < 0)//if larger than current max, update maxPos
		this.maxPosIndex = curIndex;
	}
    }

    //@dm should be --> this is the max-allowed mid-point distances for halfCluster
    public boolean withinDm(EdgeCluster eClu, int dm){
	return this.withinDm(eClu.minPosRead(), dm) && this.withinDm(eClu.maxPosRead(), dm);
    }

    
    // checking if (Read r on curEC) --> linearPosition is within DM of this edgeCluster(min & max)
    // used for leading end clustering
    // CURRENTLY NOT USED
    public boolean withinDm(Read r, int curEC, int dm){
	int toMin = r.distanceToOtherReadFromCurrentPosition(curEC
							     , readList.get(minPosIndex));
	if(toMin <= dm){
	    int toMax = r.distanceToOtherReadFromCurrentPosition(curEC
								 , readList.get(maxPosIndex));
	    if(toMax < dm)
		return true;
	}
	return false;
    }


    /* 
     * returns minimum linear poistion for this EDGE at ec
     * this will return + if the read is fwd respect to genome at ec-th interval of edge.
     * otherwise it will be a negative number.
     */
    public int getLinPosMinForGivenEC(int ec){
	if(this.minPosRead().getMappedEdge().nthInterval(ec).isFwd()){
	    if(this.minPosRead().isFwdAtEc(ec))
		return this.minPosRead().getLinPosWithEdgeCounterAt(ec) - this.minPosRead().getDistanceToEnd();
	    else
		return -1*(this.minPosRead().getLinPosWithEdgeCounterAt(ec) - this.minPosRead().getDistanceToEnd());
	}else{
	    if(this.maxPosRead().isFwdAtEc(ec))
		return this.maxPosRead().getLinPosWithEdgeCounterAt(ec) - this.maxPosRead().getDistanceToEnd();
	    else
		return -1*(this.maxPosRead().getLinPosWithEdgeCounterAt(ec) - this.maxPosRead().getDistanceToEnd());

	}
    }

    public int getLinPosMaxForGivenEC(int ec){
	if(this.maxPosRead().getMappedEdge().nthInterval(ec).isFwd()){
	    return this.maxPosRead().getLinPosWithEdgeCounterAt(ec) + this.maxPosRead().getDistanceToEnd();
	}else{
	    return this.minPosRead().getLinPosWithEdgeCounterAt(ec) + this.minPosRead().getDistanceToEnd();
	}
    }

    /* checking if two EdgeClusters are within Dm at specified edgeCoutner values. 
     * first interval's posiiton is at linPosMinFirstEdge and
     * testing interval's position(max) is computed within the method to check if the distance between these two points are within Dm
     returns -1 if false.
     returns linPosMaxRemote if true;
     */
    public int isRemoteWithinDmGivenMinPos(int linPosMinFirstEdge, boolean firstEdgeDir, int ec_remote, int dm){
	int linPosMaxRemote = 0;
	if(this.minPosRead().getMappedEdge().nthInterval(ec_remote).isFwd()){
	    if(this.maxPosRead().isFwdAtEc(ec_remote) == firstEdgeDir)
		linPosMaxRemote = this.maxPosRead().getLinPosWithEdgeCounterAt(ec_remote) + this.maxPosRead().getDistanceToEnd();
	    else
		return -1;
	}else{
	    if(this.minPosRead().isFwdAtEc(ec_remote) == firstEdgeDir)
		linPosMaxRemote = this.minPosRead().getLinPosWithEdgeCounterAt(ec_remote) + this.minPosRead().getDistanceToEnd();
	    else
		return -1;
	}
	
	if( linPosMaxRemote - linPosMinFirstEdge <= dm )
	    return linPosMaxRemote;
	return -1;
    }


    /* checking if two halfClusters are within Dm given both edgeCounters. */
    /* ec_this is alwasy smaller than ec_other */
    /*public boolean withinDm(EdgeCluster other, int ec_this, int ec_other, int dm){
        int linPosMinThis = 0;
        //int linPosMaxThis = 0;                                                                                                                                          
	//int linPosMinOther = 0;                                                                                                                                         
        int linPosMaxOther = 0;
        if(this.minPosRead().getMappedEdge().nthInterval(ec_this).isFwd()){
            linPosMinThis = this.minPosRead().getLinPosWithEdgeCounterAt(ec_this);
            //linPosMaxThis = this.maxPosRead().getLinPosWithEdgeCounterAt(ec_this);
        }else{
            linPosMinThis = this.maxPosRead().getLinPosWithEdgeCounterAt(ec_this);
            //linPosMaxThis = this.minPosRead().getLinPosWithEdgeCounterAt(ec_this);
	}
        if(other.minPosRead().getMappedEdge().nthInterval(ec_this).isFwd()){
            //linPosMinOther = other.minPosRead().getLinPosWithEdgeCounterAt(ec_other);
            linPosMaxOther = other.maxPosRead().getLinPosWithEdgeCounterAt(ec_other);
        }else{
            //linPosMinOther = other.maxPosRead().getLinPosWithEdgeCounterAt(ec_other);
            linPosMaxOther = other.minPosRead().getLinPosWithEdgeCounterAt(ec_other);
        }
	
        if( linPosMaxOther - linPosMinThis > Constatns.DM )
            return true;
        return false;
	}*/

    // checking if (Read r on with unknown EC information) is within DM of this edgeCluster(min & max)
    // CURRENTLY USED FOR :  remote end clustering
    // 
    public boolean withinDm(Read r, int dm){
	boolean toMin = Edge.isDistanceBetween2GPsWithinDM( r.getMappedEdge(), 
							    r.getGraphPos(), 
							    readList.get(minPosIndex).getMappedEdge(),
							    readList.get(minPosIndex).getGraphPos(),
							    dm);
	if(toMin){
	    boolean toMax = Edge.isDistanceBetween2GPsWithinDM( r.getMappedEdge(), 
								r.getGraphPos(), 
								readList.get(maxPosIndex).getMappedEdge(),
								readList.get(maxPosIndex).getGraphPos(),
								dm
								);
	    if(toMax)
		return true;
	}
	return false;
    }

    /*
    public GraphRange getGraphRange(){
	return new GraphRange(this.minPosRead(), this.maxPosRead());
	}*/

   
    public Read minPosRead(){
	return this.readList.get(this.minPosIndex);
    }
    
    public Read maxPosRead(){
	return this.readList.get(this.maxPosIndex);
    }
    
    public ArrayList<Read> getReadList(){
	return this.readList;
    }

    public int size(){
	return this.readList.size();
    }

    public int getMinPosIndex(){
	return this.minPosIndex;
    }

    public int getMaxPosIndex(){
	return this.maxPosIndex;
    }


    /* analogous to min-max pair */
    private int minPosIndex; //this is the index of the minPos Read in the readList
    private int maxPosIndex; //this is the index of the maxPos read in the readList
    private ArrayList<Read> readList;

    public void updateEC(int ec){
	this.ec = ec;
    }

    private int ec;

    
    //totalBases --> necessary to accurately calculate coverages.
    private int totalBases; //if readLengths is constants, it's simply numReads * readLength. But it's variable length we are dealing with, so we need to update this whenever we add read.

    

    
}
