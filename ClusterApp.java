import java.util.*;

public class ClusterApp{

    //although we have a full list of clusters, we only need to scan locally located 
    //clusters as defined by curListStart and curListEnd
    //private ArrayList<Cluster> clusterList;
    
    //this is the hashtable containing all HalfCluster Pairs. For each pair of HalfCluster,
    //only one entry is stored.
    private HashMap<Integer, HalfCluster> halfClusterTable; 

    //this is the lowerbound index boundary of clusterList for curList in algo.
    //private int curListStart; 
    
    // this is the upperbound index boundary of clusterList for curList in algo
    //private int curListEnd; 
    
    //private LinkedList<HalfCluster> activeClusterList;
    private ActiveClusterList activeClusterList; //[10/8/14]: updated as a class ActiveClusterList containing a pairing LinkedList storing integer indices for edge counter for positionReads

    private Thread t;

    //maxSize of leading/remoteEnd of a cluster ~ (insertSize-readLength)
    private int dm;// = Constants.DM - Constants.READLEN;//400; //Constants.DM 
    
    //private int MIN_NUM_READS_PER_CLUSTER = 20;

    /*

    public static void main(String[] args){
	//args[0] : threadFile
	//args[1] : sorted sam
	//args[2] : medmadfile
	//args[3] : conf file
	//ClusterApp obj = new ClusterApp(new Thread(4639675, args[0]));
	if(args.length == 4){
	    Constants.loadConstants(args[2], args[3]);
	    ClusterApp obj = new ClusterApp(new Thread(Constants.GENOMELEN, args[0]));
	    obj.run(args[1], true);
	}else{
	    System.err.println("USAGE: java CluserApp <threadFile> <sorted sam> <medmad> <config_file>");
	}
    }
    */

    public static void main(String[] args){
	//args[0] : command (depth, fullNoDepth, full)
	//args[1] : threadFile
	//args[2] : sorted sam
	//args[3] : medmadfile
	//args[4] : conf file
	//args[5] : depth data (for de-serialization)
	//ClusterApp obj = new ClusterApp(new Thread(4639675, args[0]));
	boolean error = false;
	if(args.length == 0){
	    System.err.println("USAGE: java CluserApp depth/fullNoDepth/fullMidsorted");
	}else if(args[0].equals("depth")){
	    if(args.length == 5){
		Constants.loadConstants(args[3], args[4]);
		ClusterApp obj = new ClusterApp(new Thread(Constants.GENOMELEN, args[1]));
		obj.run(args[2], true);
	    }else
		System.err.println("USAGE: java CluserApp depth <threadFile> <sorted sam> <medmad> <config_file>");
	}else if(args[0].equals("fullNoDepth")){
	    if(args.length == 6){
		Constants.loadConstants(args[3], args[4]);
		ClusterApp obj = new ClusterApp(new Thread(Constants.GENOMELEN, args[1]));
		obj.run(args[2], args[5]);
	    }else
		System.err.println("USAGE: java CluserApp fullNoDepth <threadFile> <simpleDiscordantRemoved-midpoint sorted sam> <medmad> <config_file> <depth data>");
	    
	}else if(args[0].equals("fullMidsorted")){
	    if(args.length == 5){
		Constants.loadConstants(args[3], args[4]);
		ClusterApp obj = new ClusterApp(new Thread(Constants.GENOMELEN, args[1]));
		obj.legacyRun(args[2]);
	    }else
		System.err.println("USAGE: java CluserApp fullMidSorted <threadFile> <sorted sam> <medmad> <config_file>");
	
	}else
	    System.err.println("USAGE: java CluserApp depth/fullNoDepth/fullMidsorted");
	
    }

    public void run(String samfile, boolean simple){
	System.err.println("*** Loading READS & GETTING DEPTH ARRAYS...");
	new SimpleReadsLoader(samfile, this.t);
	System.err.println("--- DONE Loading ----");
    }
    
    public void run(String discordantMidpointSortedSamFile, String depthSerializedData){
	boolean updateDepth = false;
	System.err.println("*** De-Serializing Depth Array . . .\n");
	t.deserializeDepthArray(depthSerializedData);
	System.err.println("--- DONE De-Serialization ----\n");
	
	System.err.println("*** Loading READS ...\n");
	new ReadsLoader(discordantMidpointSortedSamFile, this.t, updateDepth);
	System.err.println("--- DONE Loading ----\n");
	
	System.err.println("\n*** Removing concordant and remote-end reads ...\n");
	new ConcordantRemover(t);
	System.err.println("--- DONE removing ----\n");
	
	System.err.println("\n*** Clustering ...\n");
	this.cluster();
	this.prepHalfClustersForPairing();
	System.err.println("--- DONE Clustering ----\n");
	//this.printClusters();
	System.err.println("\n*** ASSIGN EVENTS ...\n");
	this.assignEvents();
	this.printLeftOverClusters();
	
    }

    //takes the midpoint sorted sampe file 
    //1. Loads reads that are paired (pass the quality filter)
    //2. Removes concordant read pairs
    //3. Cluster
    //4. Assign Events.
    public void legacyRun(String samfile){
	boolean updateDepth = true;
	System.err.println("*** Loading READS ...");
	new ReadsLoader(samfile, this.t, updateDepth);
	System.err.println("--- DONE Loading ----");
	System.err.println("\n*** Removing concordant and remote-end reads ...");
	new ConcordantRemover(t);
	System.err.println("--- DONE removing ----");
	System.err.println("\n*** Clustering ...");
	this.cluster();
	this.prepHalfClustersForPairing();
	//this.printClusters();
	System.err.println("\n--------- ASSIGN EVENTS ...");
	this.assignEvents();
	this.printClusters();
    }
    

    public ClusterApp(Thread t){
	this.dm = Constants.DM;//Constants.DM - 2*Constants.READLEN;
	this.t = t;
	t.traverseDemo();
	//this.clusterList = new ArrayList<Cluster>();
	this.halfClusterTable = new HashMap<Integer, HalfCluster>();
	
	//[10/8/14]
	this.activeClusterList = new ActiveClusterList();//new LinkedList<HalfCluster>(); 
	//this.curListStart = 0;
	//this.curListEnd = 0;
    }

    
    public void printClusters(){
	//StringBuffer bf = null;
	Iterator<Integer> keyItr = this.halfClusterTable.keySet().iterator();
	while(keyItr.hasNext()){
	    //bf=new StringBuffer();
	    Integer curInt = keyItr.next();
	    if(this.halfClusterTable.get(curInt).numReads() >= Constants.MIN_NUM_READS_PER_CLUSTER){
		System.err.println("##########################################");
		String thisResult = this.halfClusterTable.get(curInt).print();
		System.err.println("------------------------------------------");
		String pairResult = this.halfClusterTable.get(curInt).getPair().print();
		//if(thisResult == null || pairResult == null){
		    
		
		    //System.out.println("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&");
		    //String pairResultSorted = this.halfClusterTable.get(curInt).getPair().print(true, this.t);//true for sorting EdgeList and t is needed to sort
		System.err.println("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
		System.err.println("U>" + thisResult + "\t" + pairResult + "\n");
		    //System.out.println("S>" + thisResult + "\t" + pairResult);
		//System.out.println(bf);
		
	    }
	}
    }

    public void printLeftOverClusters(){
	System.out.println("#\t\tUNASSIGNED CLUSTERS");
	//StringBuffer bf = null;
	Iterator<Integer> keyItr = this.halfClusterTable.keySet().iterator();
	while(keyItr.hasNext()){
	    //bf=new StringBuffer();
	    Integer curInt = keyItr.next();
	    System.out.println(this.halfClusterTable.get(curInt).printAsDR());
	}
	keyItr = null;
    }


    public void cluster(){
	Edge curEdge = null;
	Interval curInterval = null;
	
	int linPos = 0;
	
	//traversing the thread --> for each interval on the thread
	for(int i =0; i < this.t.getIntervalSequence().size();i++){
	    
	    curInterval = this.t.getIntervalSequence().get(i);
	    curEdge = curInterval.getEdge();
	    
	    //gets the read and cluster list for the current edge.
	    LinkedList<Read> readList = curEdge.getReads();
	    LinkedList<HalfCluster> clusterList = curEdge.getClusters();
	    
	    //if the direction of the current interval is forward, we return the fwdIterator, otherwise reverseIterator
	    DirectionalIterator<Read> readIterator = new DirectionalIterator<Read>(readList, curInterval.isFwd());
	    DirectionalIterator<HalfCluster> clusterIterator = new DirectionalIterator<HalfCluster>(clusterList, curInterval.isFwd());
	    
	    //switches to feed read/cluster
	    boolean feedRead = true; 
	    boolean feedCluster = true;
	    
	    Read fedRead = null;
	    HalfCluster fedCluster = null;
	    
	    HalfCluster curCluster = null;
	    
	    //while we have reads or halfClusters left on the Edge.
	    while(readIterator.hasNext() || clusterIterator.hasNext()){
		if(Constants.DEBUG){
		    System.err.println("Reading reads from");
		    System.err.println(curInterval.toString() + "\n");
		}
		if(feedRead && readIterator.hasNext()){//if we need to feed nextRead
		    fedRead = readIterator.next();
		    readIterator.remove(); // since this a bare read, it is immediately removed from readList --> this will be repackaged as cluster.
		}else if(feedRead)
		    fedRead = null;
		if(feedCluster && clusterIterator.hasNext())//if we need to feed nextCluster
		    fedCluster = clusterIterator.next();
		else if(feedCluster)//if we need to feed nextCluster && no more cluster on the iterator
		    fedCluster = null;
		
		//added interval direction to reflect ordering of reverse interval.
		boolean virginCluster = this.isReadFirst( fedRead, fedCluster, curInterval.isFwd() );
		
		//if read comes first:
		if(virginCluster){
		    if(Constants.DEBUG)
			System.err.println("\tFetching a read-pair and wrapping as a cluster :\n\t\t" + fedRead.toString() + "\n\t\t" + fedRead.getPairingRead().toString());
		    curCluster = new HalfCluster(fedRead, fedRead.getPairingRead()); 
		    feedRead = true;
		    feedCluster = false;
		    //if(fedCluster == null) //newly inserted cluster needs to be recalled.
		    //	feedCluster = true;
		}else{
		    if(Constants.DEBUG)
			System.err.println("\tFetching a cluster : \n\t\t" + fedCluster.getId() + "\t[CurLinPos]:" + fedCluster.getCurLinPos());
		    curCluster = fedCluster;
		    feedCluster = true;
		    feedRead = false;
		}
		
		//get current position --> gets the position from positionRead of the cluster
		int curLinPos = this.getCurLinPos(curCluster, curInterval); 
		if(Constants.DEBUG){
		    System.err.println("\t--> curLinPos is :\t" + curLinPos + "\t(curCluLinPos:" + curCluster.getCurLinPos() + ")");
		    System.err.println("\t--> NEED TO UPDATE ActiveClusterList FIRST");
		}
		//removes any clusters farther than Dm from curLinPos from activeClusterList.
		this.updateActiveClusterList(curLinPos, curCluster);
		curCluster.setCurLinPos(curLinPos);
		if(Constants.DEBUG)
		    System.err.println("\t--> DONE UPDATING. now Setting the curLinPos. DONE SETTING :\t" + curLinPos + "\t(curCluLinPos:" + curCluster.getCurLinPos() + ")");
		
		
		
		if(this.mergeCluster(curCluster, curInterval.getEC())){
		    if(Constants.DEBUG)
			System.err.println("\tMERGED with another Cluster");
		    curCluster.resetCurLinPos();//need to set curLinPos to -1;
		    if(!virginCluster){ // since this cluster has been merged with one of the clusters on active cluster
			clusterIterator.remove(); //we remove this from the edge
			halfClusterTable.remove(new Integer(curCluster.getId())); //also remove this from the total cluster list.
		    }
		}else{//if not merged
		    //[10/8/14]
		    activeClusterList.add(curCluster, curInterval.getEC());
		    if(Constants.DEBUG)
			System.err.println("\tNOT MERGED! FIRST WE ADD curClu to the ActiveClusterList " + "\t|ActiveCluList|=" + this.activeClusterList.size());
		    if(virginCluster){ //and it's a virginCluster --> we need to add it to current Edge and the totalClusterTable
			if(Constants.DEBUG)
			    System.err.println("\t\tVirginCluster so need to ADD TO curEdge and ClusterTable.");
			if(fedCluster == null)
			    clusterIterator.add(curCluster); 
			else//this is used to insert the curCluster before the fedCluster. iterator position is already past the fedCluster so need to use this.
			    clusterIterator.addWhenNotNull(curCluster);
			
			halfClusterTable.put(new Integer(curCluster.getId()), curCluster);
			
			//activeClusterList.add(curCluster);
			
		    }
		}
	    }
	}
	
    }


    private boolean mergeCluster(HalfCluster cl, int curEC){
	boolean foundMatch = false;
	boolean foundIdentical = false;
		
	//[03/03/15] : using combined iterator to synchronize two LinkedLists used in ActiveClusterList]
	ActiveClusterListIterator iterator = this.activeClusterList.listIterator();
		
	HalfCluster clu_x = null;
	//[10/8/14]
	int clu_xEC = -1;
	if(Constants.DEBUG)
	    System.err.println("\t[MERGEcluster]:\t Attempting to merge " + cl.getId() + "\t|ActiveCluList|=" + this.activeClusterList.size());

	while(iterator.hasNext()){
	    
	    ClusterRecentECPair next = iterator.nextPair();
	    clu_x = next.getClu();
	    clu_xEC = next.getRecentEC();//[10/8/14]

	    if(Constants.DEBUG)
		System.err.println("\t\t[MERGECluster]:\t Scanning againt cluX( ID:" + clu_x.getId() +" ):\t" + cl.getId());

	    //this (running into identical cluster in activeList) should never happend but left the code as it is.
	    if(!foundIdentical 
	       && cl.getId() == clu_x.getId()){
		if(Constants.DEBUG)
		    System.err.println("\t\t[MERGECluster]:\t found identical Cluster in activeClusterList so removing the older mapping instance");
		iterator.remove(); // single remove call now remove from both iterators
		//iteratorEC.remove();
		foundIdentical = true;
		if(foundMatch)
		    break;
	    }else if(!foundMatch){
		if(Constants.DEBUG){
		    System.err.println("\t\t[MERGECluster]curClu PR:\t" + cl.getPositionRead().toString());
		    System.err.println("\t\t[MERGECluster]CLU_X PR:\t" + clu_x.getPositionRead().toString());
		}
		if(cl.hasSameLeadingOrientationAs(clu_x, curEC, clu_xEC)//[10/8/14]
		   && cl.isWithinDm(clu_x, this.dm)){
		    if(Constants.DEBUG)
			System.err.println("\t\t\t[MERGECluster]:\t curClu " + cl.getId() + " merged to clu_x:\t" + clu_x.getId());
		    
		    clu_x.addAllFrom(cl);
		    foundMatch = true;
		    //if(foundIdentical) /* commenting this out since we no longer need to check identical cluster in activeList*/
		    break;
		}
	    }
	}
	iterator = null;
	return foundMatch;
    }


    //!!!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!
    //
    //NEED TO UPDATE SO THAT THE WIDTH of CURRENT CLUSTER IS REFLECTED IN THE CALCULATION AS WELL!!!!!!
    //
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/*
    private int getCurLinPos(Cluster cl, Interval intv){
	int tmp = intv.getStart();
	if(intv.isFwd())
	    tmp += cl.getPositionRead().getPosition() -1;
	else
	    tmp += (intv.length() - cl.getPositionRead().getPosition() + 1);
	return tmp;
    }
*/
    private int getCurLinPos(HalfCluster cl, Interval intv){
	return intv.getLinearPosition(cl.getPositionRead().getGraphPos());
    }
    

    /* returns true if read at readIndex of edge comes before cluster 
     * returns false otherwise --> if tied, it returns false we want to process cluster first.
     */
    /*
    private boolean isReadFirst(int readIndex, int clusterIndex, Edge edge){
	if(readIndex >= edge.getReads().size()){
	    return false;
	}else if(clusterIndex >= edge.getReads().size()){
	    return true;
	}else{
	    int compVal = edge.getReads().get(readIndex).compareTo(edge.getClusters().get(clusterIndex).getPositionRead());
	    if(compVal < 0)
		return false;
	    else// if on same position procee cluster first
		return true;
	}
    }*/
    
    
    // Comprares curRead and curCluser to see which one comes first.
    // NOTE here curRead and curCluster are always on same Edge. 
    
    private boolean isReadFirst(Read curRead, HalfCluster curCluster, boolean intervalDirection){
	if(curRead == null)
	    return false;
	else if(curCluster == null)
	    return true;
	else{
	    //this returns compareVal based on the edge direction(direction of first interval beloning to the edge)
	    //therefore, if intervalDirection is opposite of the edge then, we return the opposite.
	    int compVal = curRead.compareToReadOnSameEdge(curCluster.getPositionRead());
	    
	    if(compVal < 0){
		if(intervalDirection)
		    return false;
		return true;
	    }else{
		if(intervalDirection)
		    return true;
		return false;
	    }
	}
    }


    // UPDATE TO REMOVE DUPLICATE clusters ( that is a copy of cluster in the activeClusterList == curClu)
    // Simply goes through the activeClusterList to remove clusters that are farther than Dm away.
    // Since the list is sorted, we immediately break once we find a cluster that is within Dm.
    private void updateActiveClusterList(int curLinPos, HalfCluster clu){
	boolean removeDuplicate = true;
	boolean foundBoundary = false;
	if(clu.getCurLinPos() == -1) // curLinPos is -1. we know this cluster is not in the activeClusterList.
	    removeDuplicate = false;
	
	int minPos = curLinPos - dm + 1; // this is the new cutff which is within dm away from current position
	if(Constants.DEBUG)
	    System.err.println("\t\t@@@ UPDATING ActiveClusterList with minPos = " + minPos);
	
	ActiveClusterListIterator iterator = this.activeClusterList.listIterator(true);
	HalfCluster curClu = null;
	int x = 0;
	while(iterator.hasNext()){
	    curClu = iterator.next();
	    if(Constants.DEBUG)
		System.err.println("\t\t@@@ checking CLU_X( " + curClu.getId() + "(" + x + ") ");
	    if(foundBoundary && !removeDuplicate){//we break as soon as we are done removing duplciate && finding boundary
		if(Constants.DEBUG)
		    System.err.println("\t\t\tWe already found duplicates and boundary: so skipping");
		break;
	    }
	    else if(removeDuplicate && (clu.getId() == curClu.getId())){//if we haven't found duplicate, remove if it's a duplicate
		if(Constants.DEBUG)
		    System.err.println("\t\t\tFOUND DUPLICATE. REMOVING IT from activeClusterList");
		curClu.resetCurLinPos(); // this is unecessary as curLinPos will be updated after the call to this method. but just keeping it consistent
		if(Constants.DEBUG)
		    System.err.println("\t\t\t\tHas the linPos rest properly? : " + curClu.getCurLinPos() + ":" + clu.getCurLinPos());
		iterator.remove(); 
		removeDuplicate = false;
	    }else if(curClu.isWithinBoundary(minPos)){
		if(Constants.DEBUG)
		    System.err.println("\t\t\tWITHIN BOUNDARY SO NO NEED TO FURTHER CHECK other than duplicateRemoval");
		foundBoundary = true;
	    }else{
		if(Constants.DEBUG)
		    System.err.println("\t\t\tNOT WITHIN BOUNDARY SO REMOVING IT");
		curClu.resetCurLinPos();
		if(Constants.DEBUG)
		    System.err.println("\t\t\t\tHas the linPos rest properly? : " + curClu.getCurLinPos() + ":" + clu.getCurLinPos());
		iterator.remove();
	    }
	    x++;
	}
	iterator = null;
    }

    
    /* binrary search to find the lowerbound on the 5' side of the leading end of cluster */
    /* special cases (2) are handled prior to binary search */
    /*
    private void updateCurListBoundary(int curLinPos){
		
	int minPos = curLinPos - dm + 1; // this is the new cutoff which is within dm away from current position
	
	if(clusterList.get(curListStart).getLinearPosition() >= minPos)
	    this.curListEnd++; 
	else if(clusterList.get(curListEnd).getLinearPosition() < minPos)
	    this.curListStart = curListEnd + 1;
	else{
		   
	    int matchIndex = (curListStart + curListEnd) / 2;
	    while(clusterList.get(matchIndex).getLinearPosition() != minPos
		  && curListStart <= matchIndex){
		if(clusterList.get(matchIndex).getLinearPosition() < minPos)
		    matchIndex = (matchIndex + curListEnd) /2;
		else
		    matchIndex = (curListStart + matchIndex) /2;
	    }
	    
	    this.curListStart = matchIndex; // start position is updated.
	    this.curListEnd++; // upper bound needs to be increased by 1 since this is index which now points to the current cluster
	}
    }*/



    
    private void prepHalfClustersForPairing(){
	/*	Iterator<Integer> keyItr = this.halfClusterTable.keySet().iterator();
	//ArrayList<Integer> 
	while(keyItr.hasNext()){
	    Integer curInt = keyItr.next();
	    HalfCluster tmp = this.halfClusterTable.get(curInt);
	    if(tmp.numReads() >= Constants.MIN_NUM_READS_PER_CLUSTER){
		System.err.println("Computing DirectionalRangesPair for :\t" + tmp.getId() + "\t" + tmp.getPair().getId());
		tmp.computeDirectionalRangesPair();
	    }//else
	    //	this.halfClusterTable.remove(curInt);
	    }*/

	Iterator<Map.Entry<Integer,HalfCluster>> itr = this.halfClusterTable.entrySet().iterator();
	while(itr.hasNext()){
	    Map.Entry<Integer,HalfCluster> entry = itr.next();
	    HalfCluster tmp = entry.getValue();//this.halfClusterTable.get(curInt);
	    if(tmp.numReads() >= Constants.MIN_NUM_READS_PER_CLUSTER){
		//System.err.println("Computing DirectionalRangesPair for :\t" + tmp.getId() + "\t" + tmp.getPair().getId());
		tmp.computeDirectionalRangesPair();
	    }else
	    	itr.remove();//this.halfClusterTable.remove(curInt);
	}
    }
    

    

    private void assignEvents(){
	//this.prepHalfClustersForPairing();

	/* 0. get list of depth depleted ranges */
	ArrayList<DirectionalRange> deletionRangeList = this.t.findDepthDepletedRegions();
	

	/* 1. screen for simple deletion or deletion via homologus recombination*/
	HashMap<Integer, Integer> assignedHalfClusterIDs = new HashMap<Integer, Integer>();

	StringBuffer deletionBuffer = new StringBuffer();
	StringBuffer transpositionBuffer = new StringBuffer();
	StringBuffer inversionBuffer = new StringBuffer();
	StringBuffer tandemBuffer = new StringBuffer();
	
	Iterator<Map.Entry<Integer,HalfCluster>> itr = this.halfClusterTable.entrySet().iterator();
	while(itr.hasNext()){
	    Map.Entry<Integer,HalfCluster> entry = itr.next();
	    HalfCluster hcl1 = entry.getValue();//this.halfClusterTable.get(curInt);
	    boolean debug = false;
	    /*if(hcl1.getId() == 8493)
	      debug = true;*/
	    //returns a length-3 int array. 
	    //int[0] : 0 if  ---------this----->       <----this'spair--- deletion
	    //       : 1 if ----this's pair---->       <------this-----   deletion
	    //       : -1 if no deletion.
	    //int[1] : index for accesing directionalRange for 5'side i
	    //int[2] : index for accesing directionalRange for 3'side i
	    if(debug)
		System.err.println("======== CHECKING DELETION FOR 14631");
	    int[] pairIndicies = hcl1.isSimpleDeletion(this.t);
	    if(debug){
		System.err.println("======== 14631: pairIndicies[0] : " + pairIndicies[0]);
		System.err.println("======== hcl1 checkAgain Flag\t" + hcl1.checkRepDeletionAgain() );
	    }
	    if(pairIndicies[0] > -1){
		deletionBuffer.append(hcl1.getDeletionEventString(pairIndicies));
		if(Constants.DEBUG)
		    hcl1.printDirect();
		itr.remove();//this.halfClusterTable.remove(keys[i]);
		hcl1.removeDepletedRegionsThatMatchDeletionSign(deletionRangeList, pairIndicies);
		assignedHalfClusterIDs.put(new Integer(hcl1.getId()), new Integer(hcl1.getId()));
	    }//else if(pairIndicies[0] == -1)
		
	}

	itr = null;


	/* 2. screen for TRANSPOSITION & INVERSIONS */
	Integer[] keys = this.halfClusterTable.keySet().toArray(new Integer[1]);
	ArrayList<Integer> list2remove = new ArrayList<Integer>();


	boolean debug = false;
	
	for(int i=0;i<keys.length;i++){
	    HalfCluster hcl1 = this.halfClusterTable.get(keys[i]);
	    if(assignedHalfClusterIDs.get(new Integer(hcl1.getId())) == null){
		for(int j=(i+1); j< keys.length;j++){
		    HalfCluster hcl2 = this.halfClusterTable.get(keys[j]);
		    
		    if(assignedHalfClusterIDs.get(new Integer(hcl2.getId())) == null){
			//if( (hcl1.getId() == 18509 && hcl2.getId() == 18551) || (hcl2.getId() == 18551 && hcl1.getId() == 18509))
			//    debug = true;
			//else
			//    debug = false;
			String eventType = hcl1.getEventType(hcl2, debug);
			if(eventType != null){
			    int bufferdId = ClusterApp.getEventTypeID(eventType);
			    if(Constants.DEBUG){
				System.err.println("****HCL1***");
				hcl1.printDirect();
				System.err.println("****HCL2***");
				hcl2.printDirect();
			    }


			    list2remove.add(new Integer(hcl1.getId()));
			    list2remove.add(new Integer(hcl2.getId()));
			    assignedHalfClusterIDs.put(new Integer(hcl1.getId()), new Integer(hcl2.getId()));
			    assignedHalfClusterIDs.put(new Integer(hcl2.getId()), new Integer(hcl1.getId()));
			    break;
			    /*if(bufferdId == 0 || bufferId == 1)
			      transpoistionBuffer.append(hcl1.getEventString(hcl2, eventType));
			      else if(bufferId = 2)
			      inversionBuffer.append(hcl1.getEventString(hcl2, eventType));*/
			}
		    }
		}
	    }
	}
	
	//removes the leading cluster from the clusterTable.
	for(int i=0; i<list2remove.size(); i++)
	    this.halfClusterTable.remove(list2remove.get(i));
	list2remove = null;
	
	
	
	
	/*REPETATIVE DELETION & TANDEM DUPLICATION HERE*/
	itr = this.halfClusterTable.entrySet().iterator();
	while(itr.hasNext()){
	    Map.Entry<Integer,HalfCluster> entry = itr.next();
	    HalfCluster hcl1 = entry.getValue();//this.halfClusterTable.get(curInt);
	    //System.err.println("hcl1.Id:\t" + hcl1.getId());
	    debug = false;
	    /*if(hcl1.getId() == 8493){
		System.err.println("14631 FOUND");
		debug = true;
	    }
	    if(hcl1.getId() == 8494){
		System.err.println("14632 FOUND");
		debug = true;
		}*/
	    if(hcl1.checkRepDeletionAgain()){
		//System.err.println("\tChecking repetative deletion");
		//returns a length-3 int array. 
		//int[0] : 0 if  ---------this----->       <----this'spair--- deletion
		//       : 1 if ----this's pair---->       <------this-----   deletion
		//       : -1 if no deletion.
		//int[1] : index for accesing directionalRange for 5'side i
		//int[2] : index for accesing directionalRange for 3'side i
		int[] pairIndicies = hcl1.isSimpleDeletion(this.t);
		if(debug)
		    System.err.println("pairIndicies[0] : \t" + pairIndicies[0]);
	       
		if(pairIndicies[0] == -1 || pairIndicies[0] == -2){
		    deletionBuffer.append(hcl1.getDeletionEventString(pairIndicies) + "\n");
		    if(Constants.DEBUG)
			hcl1.printDirect();
		    itr.remove();//this.halfClusterTable.remove(keys[i]);
		    hcl1.removeDepletedRegionsThatMatchDeletionSign(deletionRangeList, pairIndicies);
		    assignedHalfClusterIDs.put(new Integer(hcl1.getId()), new Integer(hcl1.getId()));
		}
	    }else{
		int[] pairIndicies = hcl1.isTandemDuplication();
		if(pairIndicies[0] == 0 || pairIndicies[0] == 1){
		    tandemBuffer.append(hcl1.getTandemDuplicationEventString(pairIndicies) + "\n");
		    itr.remove();
		    assignedHalfClusterIDs.put(new Integer(hcl1.getId()), new Integer(hcl1.getId()));
		}else if(pairIndicies[0] == -2){
		    itr.remove();
		}
	    }
	}
	
	/* checking for recombination deletion of 2 existing copies (both flanks:repeat-bounded*/
	for(int i=0; i<deletionRangeList.size(); i++){
	    int[] updatedRange = this.checkIfFlankersAreRepeat(deletionRangeList.get(i));
	    if(updatedRange == null)
		System.out.println(deletionRangeList.get(i).toDeletionRange());
	    else
		System.out.println(new DirectionalRange(updatedRange[0], updatedRange[1]).toDeletionRangeRepeat(updatedRange));
	}
	
	
	//System.err.println("==============LEFT OVERS============");

	//for(int i=0;i<keys.length;i++){
	//    HalfCluster hcl = this.halfClusterTable.get(keys[i]);
	//    hcl.printPlain();
	//}
    }

    
    //returns null if not repeat bounded
    //returns int[6] 
    //        deletion boundary  (int[0], int[1])
    //        5' repeat boundary (int[2], int[3])
    //        3' repeat boundary (int[4], int[5])
    public int[] checkIfFlankersAreRepeat(DirectionalRange dr){
	boolean debug = false;
	if(dr.getMin() == 2065377){
	    debug = true;
	}

	int maxRepeatLength = 5000;

	//    st****min--------DELETED-----------max****end
	int min = dr.getMin() - 50;
	int max = dr.getMax() + 50;
	
	int st = min - maxRepeatLength;
	
	int end = max + maxRepeatLength;

	int[] maxFallPositions = t.doesNOverlapOnRange(max, st, min);
	
	int[] minFallPositions = t.doesNOverlapOnRange(min, max, end);

	int rpst1 = 0;
	int rpend1 = 0;
	int rpst2 = 0;
	int rpend2 = 0;
	
	if(debug){
	    System.err.println("min:\t" + min + "(" + dr.getMin() + ")");
	    System.err.println("min:\t" + max + "(" + dr.getMax() + ")");
	    System.err.println("st:\t" + st);
	    System.err.println("end:\t" + end);
	    for(int i=0; i<maxFallPositions.length;i++){
		if(i==maxFallPositions.length-1)
		    System.err.print("\t");
		System.err.println("maxFallPositions[" + i + "]:\t" + maxFallPositions[i]);
	    }
	    for(int i=0; i<minFallPositions.length;i++){
		if(i==0)
		    System.err.print("\t");
		System.err.println("minFallPositions[" + i + "]:\t" + minFallPositions[i]);
	    }
	}

	boolean found = false;
	
	if(maxFallPositions.length > 0 
	   && minFallPositions.length > 0){
	    int repeat5Len = maxRepeatLength - maxFallPositions[maxFallPositions.length - 1];
	    int repeat3Len = minFallPositions[0];
 	    int diff = repeat5Len - repeat3Len;
	    if(debug){
		System.err.println("repeat5Len:\t" + repeat5Len);
		System.err.println("repeat3Len:\t" + repeat3Len);
		System.err.println("diff:\t" + diff);
	    }
	    
	    if(repeat5Len > 50 && repeat3Len > 50 && diff >= -90 && diff <=90){
		if(debug){
		    System.err.println("#PASS the FIRST CONDITION: repeat5Len > 50 && repeat3Len > 50 && diff >= -90 && diff <=90");
		}
		//found = true;
		int[] threadIndicies5prime = this.t.getIntervalsForRange(st + maxFallPositions[maxFallPositions.length - 1] - 1, min);
		int[] threadIndicies3prime = this.t.getIntervalsForRange(max, max + minFallPositions[0] - 1);
		if(debug){
		    for(int i=0;i<threadIndicies5prime.length; i++)
			System.err.print("TI5["+ i + "]=" + threadIndicies5prime[i] + "\t");
		    System.err.println();
		    for(int i=0;i<threadIndicies3prime.length; i++)
			System.err.print("TI3["+ i + "]=" + threadIndicies3prime[i] + "\t");
		    System.err.println();
		}
		if(threadIndicies5prime.length == threadIndicies3prime.length){
		    if(debug)
			System.err.println("#SAME NUMBER OF INTERVALS THROUGH PATH");
		    for(int i=0; i<threadIndicies5prime.length; i++){
			if(t.getNthInterval(threadIndicies5prime[i]).isSameDirectionIntervalOfSameEdge(t.getNthInterval(threadIndicies3prime[i])))
			    found = true;
		    }
		}
		if(found){
		    
		    IntervalPos minIntervalPos = t.getIntervalPosForPositionX(min);
		    int distanceToEnd5 = minIntervalPos.getInterval().distanceToNextIntervalFromLinPos(min);
		    
		    IntervalPos minFallIntervalPos = t.getIntervalPosForPositionX(max + minFallPositions[0] -1);
		    int distanceToEnd3 = minFallIntervalPos.getInterval().distanceToNextIntervalFromLinPos(max + minFallPositions[0] - 1);
		    
		    if(debug){
			System.err.println("#EXTENSION");
			System.err.println("distanceToEnd5:\t" + distanceToEnd5);
			System.err.println("distanceToEnd3:\t" + distanceToEnd3);
		    }
		    int jump = 1;
		    while(true){
			if( t.getNthInterval(minIntervalPos.getThreadIndex()+jump).isSameDirectionIntervalOfSameEdge(t.getNthInterval(minFallIntervalPos.getThreadIndex()+jump)) ){
			    if(debug)
				System.err.println("NEXT INTERVALS MATCH! EXTENSION:\t" + jump );
			    distanceToEnd5 += t.getNthInterval(minIntervalPos.getThreadIndex()+jump).length();
			    distanceToEnd3 += t.getNthInterval(minFallIntervalPos.getThreadIndex()+jump).length();
			    jump++;
			}else
			    break;
			
		    }

		    IntervalPos maxFallIntervalPos = t.getIntervalPosForPositionX(st + maxFallPositions[maxFallPositions.length - 1] -1);
		    int distanceToSt5 = maxFallIntervalPos.getInterval().distanceToPreviousIntervalFromLinPos(st + maxFallPositions[maxFallPositions.length - 1] -1);

		    IntervalPos maxIntervalPos = t.getIntervalPosForPositionX(max);
		    int distanceToSt3 = maxIntervalPos.getInterval().distanceToPreviousIntervalFromLinPos(max);

		    if(debug){
			System.err.println("#EXTENSION");
			System.err.println("distanceToSt5:\t" + distanceToSt5);
			System.err.println("distanceToSt3:\t" + distanceToSt3);
		    }

		    jump = 1;
		    while(true){
			if( t.getNthInterval(minIntervalPos.getThreadIndex()-jump).isSameDirectionIntervalOfSameEdge(t.getNthInterval(minFallIntervalPos.getThreadIndex()-jump)) ){
			    if(debug)
				System.err.println("NEXT INTERVALS MATCH! EXTENSION:\t" + "(-)"+jump );
			    distanceToSt5 += t.getNthInterval(minIntervalPos.getThreadIndex()-jump).length();
			    distanceToSt3 += t.getNthInterval(minFallIntervalPos.getThreadIndex()-jump).length();
			    jump++;
			}else
			    break;
			
		    }
		    
		    rpst1 = st + maxFallPositions[0] - 1 - distanceToSt5;
		    min = min + distanceToEnd5;
		    rpend1 = min;

		    rpst2 = max - distanceToSt3;
		    max = max + minFallPositions[0] - 1 + distanceToEnd3;
		    rpend2 = max;
		}
	    }
	}

	if(!found)
	    return null;
	else{
	    int[] nums = new int[6];
	    nums[0] = min;
	    nums[1] = max;
	    nums[2] = rpst1;
	    nums[3] = rpend1;
	    nums[4] = rpst2;
	    nums[5] = rpend2;
	    return nums;
	}
    }
    /*
    public boolean checkIfFlankersAreRepeat(DirectionRange dr){
	int checkUpToDistance = 5000;
	int trimSize= 100;
	
	int min = dr.getMin() - trimSize;
	int max = dr.getMax() - trimSize;
	
	// 5' side 
	int start = min - checkUptoDistance;
	
	int startThreadIndex = t.getIntervalForPositionX(start).getThreadIndex();
	int minThreadIndex = t.getIntervalForPositionX(min).getThreadIndex();
	
	// 3' side 
	int end = max + checkUptoDistance;
	
	int maxThreadIndex = t.getIntervalForPositionX(max).getThreadIndex();
	int endThreadIndex = t.getIntervalForPositionX(end).getThreadIndex();
	
	ArrayList<Integer> matchingPoints = new ArrayList<Integer>();
	for(int i=start; i<=min; i++){
	    if(t.getNthInterval(i).isSameDirectionIntervalOfSameEdge(t.getNthInterval(maxThreadIndex)))
		matchingPoints.add(new Integer(i));
	}
	if(matchPoints.size() > 0){
	    
	}else
	    return false;
	
	    }*/

    /*
    public boolean checkIfFlankersAreRepeat(DirectionalRange dr){
	int checkUptoDistance = 5000;
	int min = dr.getMin();
	int max = dr.getMax();
	
	// 5' side 
	int start = min - checkUptoDistance;
	
	int startThreadIndex = t.getIntervalForPositionX(start).getThreadIndex();
	int minThreadIndex = t.getIntervalForPositionX(min).getThreadIndex();
	
	// 3' side 
	int end = max + checkUptoDistance;
	
	int maxThreadIndex = t.getIntervalForPositionX(max).getThreadIndex();
	int endThreadIndex = t.getIntervalForPositionX(end).getThreadIndex();

	
	if(!t.getNthInterval(maxThreadIndex).sufficientDistanceToNextInterval(max)){
	    maxThreadIndex++;
	    max = t.getNthInterval(maxThreadIndex).getStart();
	}
	if(!t.getNthInterval(minThreadIndex).sufficientDistanceToPrevInterval(min)){
	    minThreadIndex--;
	    min = t.getNthInterval(minThreadIndex).getStart() + t.getNthInterval(minThreadIndex).length() - 1;
	}
	
	boolean flag = false;
	//sweeping the 5'side from the start
	for(int i=startThreadIndex; i<=minThreadIndex; i++){
	    flag = false;
	    if(t.getNthInterval(i).isSameDirectionIntervalOfSameEdge(t.getNthInterval(maxThreadIndex))){
		//flag = false;
		if(i==startThreadIndex){
		    if( (t.getNthInterval(maxThreadIndex).linPos2IntervalPos(max) + 100) >= t.getNthInterval(i).linPos2IntervalPos(start) )
			flag = true;
		}else{
		    if( (t.getNthInterval(maxThreadIndex).linPos2IntervalPos(max)  + 100) >= t.getNthInterval(i).linPos2IntervalPos(t.getNthInterval(i).getStart()) )
			flag = true;
		}
		if(flag){
		    int nextI = i+1;
		    if(nextI <= minThreadIndex){
			for(int j=maxThreadIndex+1; j<= endThreadIndex; j++){
			    if(t.getNthInterval(nextI).isSameDirectionIntervalOfSameEdge(t.getNthInterval(j))){
				if(nextI == minThreadIndex){
				    flag = false;
				    if(j==endThreadIndex){
					if( t.getNthInterval(minThreadIndex).linPos2IntervalPos(min) <= (t.getNthInterval(j).linPos2IntervalPos(end) + 100) )
					    return true;
				    }else{
					if( t.getNthInterval(minThreadIndex).linPos2IntervalPos(min) <= (t.getNthInterval(j).linPos2IntervalPos(t.getNthInterval(j).getEnd()) + 100) )
					    return true;
				    }
				}else
				    nextI++;
			    }else
				break;
			}
		    }
		}
	    }
	}
	return false;
	}*/

    //return 0 : directed transposition
    //return 1 : inverted transposition
    //return 2 : inversion
    //return -1 : otherwise
    public static int getEventTypeID(String eventStr){
	if(eventStr.startsWith("1-") || eventStr.startsWith("2-") || eventStr.startsWith("3-") || eventStr.startsWith("4-"))
	    return 0;
	else if(eventStr.startsWith("5-1") || eventStr.startsWith("5-2") || eventStr.startsWith("5-3") || eventStr.startsWith("5-4") ||
		eventStr.startsWith("6-1") || eventStr.startsWith("6-2") || eventStr.startsWith("6-3") || eventStr.startsWith("6-4") )
	    return 1;
	else if(eventStr.startsWith("5-A") || eventStr.startsWith("5-B") || eventStr.startsWith("5-C") || eventStr.startsWith("5-D") ||
		eventStr.startsWith("6-A") || eventStr.startsWith("6-B") || eventStr.startsWith("6-C") || eventStr.startsWith("6-D") )
	    return 2;
	else if(eventStr.startsWith("I-1") || eventStr.startsWith("I-2") || eventStr.startsWith("I-3") || eventStr.startsWith("I-4") ||
		eventStr.startsWith("I-5") || eventStr.startsWith("I-6") || eventStr.startsWith("I-7") || eventStr.startsWith("I-8") )
	    return 2;
	else
	    return -1;
    }
    
    /*
    private void asignEvents(){
	
	//Value is:                   -1   -->  for single signal event such simple deletion or deletion via homologuous recombination
	//          pairingHalfClusterID  -->  for double signal event.
	HashMap<Integer, Integer> assignedHalfClusterIDs = new HashMap<Integer, Integer>();
	
	// 1. screen duplicative transposition 
	//Iterator<Integer> keyItr = this.halfClusterTable.keySet().iterator();
	Integer[] keys = this.halfClusterTable.keySet().toArray(new Integer[1]);
	for(int i=0;i< keys.length;i++){
	    HalfCluster hcl1 = this.halfClusterTable.get(keys[i]);
	    for(int j=(i+1); j< keys.length;j++){
		HalfCluster hcl2 = this.halfClusterTable.get(keys[j]);
		ArrayList<int[]> indexPairList_1_2_withinFacing = hcl1.isWithinFacing(hcl2);
		if(indexPairList_1_2_withinFacing != null){
		    ArrayList<int[]> ipList_1pair_2pair_withinFacing = hcl1.getPair().isWithinFacing(hcl2.getPair());
		    //checking for inversion
		    if(indexPairList_1_2_withinFacing != null){
			System.out.println("INVERSION:\t" + hcl1.print() + "\t::::\t" + hcl2.print());
			assignedHalfClusterIDs.put(keys[i], keys[j]);
			//assignedHalfCluster
		    }
		    
		    ArrayList<int[]> ipList_1_2pair_boundaryDefiner = hcl1.getPair().isBoundaryDefiner(hcl2.getPair());
		      if(indexPairList_1_2pair_boundaryDefiner != null){
		      System.out.println("TRANSPOSITION:\t" + hcl1.print() + "\t::::\t" + hcl2.print());
		      }
		}
	    }
	}
    }*/

}
