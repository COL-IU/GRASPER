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

public class ConcordantRemover{

    private HashMap<String, AlignList> alignHash;

    public ConcordantRemover(Thread t){
	this.alignHash = new HashMap<String, AlignList>();
	this.removeConcordantPairs(t);
    }
    
    private void removeConcordantPairs(Thread t){
	Interval curInterval = null;
	Edge curEdge = null;
	Read curRead = null;
	AlignList curAlignList = null;
	
	for(int i=0; i<t.getIntervalSequence().size();i++){
	    curInterval = t.getIntervalSequence().get(i);
	    curEdge = curInterval.getEdge();
	    
	    LinkedList<Read> readList = curEdge.getReads();

	    DirectionalIterator<Read> readIterator = new DirectionalIterator<Read>(readList,curInterval.isFwd());
	    
	    while(readIterator.hasNext()){
		
		curRead = readIterator.next();
		if(Constants.DEBUG2){
		    System.err.println("Accessing reads on interval[" + curInterval.getEC() + "] of edge " + curEdge.getId() + " | Lengths[inerval|edge] = [" + curInterval.length() + "|" + curEdge.length() + "]");
		    System.err.println("TOPLOOP:\t" + curRead.toStringAtEc(curInterval.getEC()));
		}
		//read was marked to be deleted, delete
		if(curRead.delete()){
		    if(Constants.DEBUG2)
			System.err.println("\tremove(1): read marked <DELETE>");
		    //curRead.demarkDelete();
		    readIterator.remove();
		}
		//if fwd
		else if(curRead.isFwdAtEc(curInterval.getEC())){
		    if(Constants.DEBUG2)
			System.err.println("\tRead curDirecton is fwd so checking alignList");
		    //if alignHash exists update the list and add the current mapping instance
		    if( (curAlignList = alignHash.get(curRead.getReadName())) != null ){
			if(Constants.DEBUG2)
			    System.err.println("\t\tAlignList exists! Updating the alignList");
			curAlignList.updateAlignList(curRead, curInterval.getEC());
		    }else{ //otherwise create new AlignList
			if(Constants.DEBUG2)
			    System.err.println("\t\tNO AlignList. Making a fresh alignList.");
			alignHash.put(curRead.getReadName(), new AlignList(curRead, curInterval.getEC()));
		    }
		}
		//if rev
		else{ 
		    if(Constants.DEBUG2)
			System.err.println("\tRead curDirection is rev so checking concordancy.");
		    if( (curAlignList = alignHash.get(curRead.getReadName())) != null){
			
			//if there are any concordant mapping instance pairing with the reverse read
			//if concordant, this call REMOVES the mapping instance from alignList
			if(curAlignList.isConcordant(curRead, curInterval.getEC())){
			    try{
				alignHash.remove(curRead.getReadName());//remove the read from hash.
				curRead.getPairingRead().markDelete();//mark the pair to be deleted.
				if(Constants.DEBUG2)
				    System.err.println("\t\tConcordant means we no longer need to keep alignHash for curRead!");
				readIterator.remove();//remove the read from the edge.
			    }catch(Exception e){
				System.err.println(curRead.toString());
			    }
			}
		    }
		}//end of else
	    }//end of while
	    readIterator = null;
	}//end of for
	
	this.alignHash = null;
	curInterval = null;
	curEdge = null;
	curRead = null;
	curAlignList = null;
	
	//
	// AT THIS POINT, These types of reads are on graph
	// 1) Discordant read pairs (both reads) --> not makred 'DELETE'
	// 2) Concordant read pairs (MAX OF 1 READ IS ON THE GRAPH --> should be marked 'DELETE'
	// 3) Singleton read
	//
	// REMOVES LEFT OVER concordant read pairs
	// AFTER this, we are left just with discordant read pairs
	// 
	// THIS ALSO removes pairing read of all discordant read-pairs
	for(int i=0; i<t.getIntervalSequence().size();i++){
	    curInterval = t.getIntervalSequence().get(i);
	    curEdge = curInterval.getEdge();
	    
	    LinkedList<Read> readList = curEdge.getReads();
	    DirectionalIterator<Read> readIterator = new DirectionalIterator<Read>(readList,curInterval.isFwd());
	    
	    while(readIterator.hasNext()){
		curRead = readIterator.next();
		if(curRead.delete()){ //if curRead is marked 'DELETE'
		    if(curInterval.isLastInterval()){
			//curRead.demarkDelete(); //no need to demark.
			if(Constants.DEBUG2){
			    System.err.println("remove(3)\t" + curRead.toString());
			}
			readIterator.remove();
		    }
		}else{//if curRead is not marked deleted, we mark the pair deleted right away to ensure deleting of remote pair.
		    //before doing so we check if a read has pair or not.
		    //if it doesn't have a pair we need to remove it at the end of the traversal path.
		    
		    if(curRead.getPairingRead() == null){//IF SINGLETON
			if(curInterval.isLastInterval()){
			    curRead.markDelete();
			    if(Constants.DEBUG2)
				System.err.println("remove(4):\t" + curRead.toString());
			    readIterator.remove();
			}
		    }
		    //IF NOT SINGLETON, THEN IT MUST BE DISCORDANT READ
		    //IN THE CASE OF DISCORDANT READ, ONLY WANT TO REMOVE whichever read that comes later.
		    //Here, we check whether the pairingRead has been marked 'DELETE' or not. IF NOT MARKED DELETE (BOTH NOT MARKED DELETE) 
		    //curRead must be the read that comes first in the thread than its pairing read.
		    else if(!curRead.getPairingRead().delete()){
			if(Constants.DEBUG2)
			    System.err.println("we are marking the pair to be deleted:\t" + curRead.getPairingRead().toString());
			curRead.getPairingRead().markDelete();
		    }
		}
	    }
	}
	
	
    }
    
}

class AlignList{
    
    //private Read r;
    private String rName;
    private LinkedList<Integer> alignList;
    private LinkedList<Boolean> boolList; //adding to distinguish read 1 and read 2

    AlignList(){
	this.rName = null;
	this.alignList = new LinkedList<Integer>();
	this.boolList = new LinkedList<Boolean>();
    }


    AlignList(Read rd, int ec){
	this();
	this.rName = rd.getReadName();
	this.addMappingInstance(rd, rd.getLinPosWithEdgeCounterAt(ec));
	/*
	if(rd.isFwdAtEc(ec)){
	    this.r = rd;
	    this.alignList = new LinkedList<Integer>();
	    this.alignList.add(new Integer(this.r.getLinPosWithEdgeCounterAt(ec)));
	    
	    this.boolList = new LinkedList<Boolean>();
	    this.boolList.add( (this.r.isFirstRead() ? Boolean.TRUE : Boolean.FALSE) );
	    }*/
    }

    private void addMappingInstance(Read rd, int linPos){
	this.alignList.add(new Integer(linPos));
	this.boolList.add( (rd.isFirstRead() ? Boolean.TRUE : Boolean.FALSE) );
    }

    public void updateAlignList(Read rd, int ec){
	int curLinPos = rd.getLinPosWithEdgeCounterAt(ec);
	ListIterator<Integer> itr = this.alignList.listIterator(0);
	ListIterator<Boolean> bitr = this.boolList.listIterator(0);
	while(itr.hasNext()){
	    int tmp = itr.next().intValue();
	    boolean read1 = bitr.next().booleanValue();
	    if(Constants.DEBUG2){
		System.err.println("\t\t\t[updateAlignList] distance : [" + curLinPos + " - " + tmp + "] = " + (curLinPos - tmp));
		System.err.println("\t\t\t[updateAlignList] MAX_D_MIDPOINTS_CORDANT(rd) : " + Constants.compute_MAX_D_MIDPOINTS_CORDANT(rd));
	    }
	    if(curLinPos - tmp > Constants.compute_MAX_D_MIDPOINTS_CORDANT(rd)){
		if(Constants.DEBUG2)
		    System.err.println("\t\t\t[updateAlignList]REMOVE: farther than max distance allowed!");
		itr.remove();
		bitr.remove();
	    }else{
		if(Constants.DEBUG2)
		    System.err.println("\t\t\t[updateAlignList]Within: so breaking.");
		break; //list is sorted so we immediately break as soon as we see any read that is within MAX_D
	    }
	}
	this.addMappingInstance(rd, curLinPos);//add new mapping instance.
    }
    
    public boolean isConcordant(Read revRd, int ec){
	int curLinPos = revRd.getLinPosWithEdgeCounterAt(ec);
	ListIterator<Integer> itr = this.alignList.listIterator(0);
	ListIterator<Boolean> bitr = this.boolList.listIterator(0);
	while(itr.hasNext()){
	    int tmp = itr.next().intValue();
	    boolean read1 = bitr.next().booleanValue();
	    
	    if(read1 != revRd.isFirstRead()){//we only pair read1 and read2. Can't pair read1-read1, read2-read2
		int dist = curLinPos - tmp;
		if(Constants.DEBUG2){
		    System.err.println("\t\t\t[isConcordant] distance : " + dist);
		}
		//since alignList is sorted, as we iterater the value of dist gets smaller
		//so once dist is smaller than MIN_D, we return false
		if(dist < Constants.compute_MIN_D_MIDPOINTS_CORDANT(revRd)){
		    if(Constants.DEBUG2){
			System.err.println("\t\t\t\t[isConcordant] distance is shorter than MIN[" + Constants.compute_MIN_D_MIDPOINTS_CORDANT(revRd) + "] so breaking.");
		    }
		    return false;
		}else if(dist <= Constants.compute_MAX_D_MIDPOINTS_CORDANT(revRd)){//else if(dist <= Constants.MAX_D_MIDPOINTS){
		    if(Constants.DEBUG2){
			System.err.println("\t\t\t\t[isConcordant]REMOVE : concordant! so removing! distance is shorter than MAX[" +Constants.compute_MAX_D_MIDPOINTS_CORDANT(revRd) + "]");
		    }
		    itr.remove();
		    bitr.remove();
		    return true;
		}
	    }
	}
	return false;
    }
    
}


