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

//indicies for each EDGE in edgeList
public class ECIndicies{

    private int len;

    /* length for these arrays are same = #Edges in edgeList*/
    private int[] curIndicies; //this stores the current threadIndexForIntervals for each edge.
    private int[] edgeIndicies; //this stores the index to access the corresponding Edge in edgeList of HalfCluster. 
    private int[] nextECs; //this stores edgeCounters for accesing next interval from Edge
    
    /* tmp indicies needed when sorting */
    private int[] curIndiciesTmp;
    private int[] edgeIndiciesTmp;
    private int[] nextECsTmp;

    public ECIndicies(ArrayList<Edge> edgeList){
	
	this.len = edgeList.size();

	this.curIndicies = new int[len];
	this.edgeIndicies = new int[len];
	this.nextECs = new int[len];

	this.curIndiciesTmp = new int[len];
	this.edgeIndiciesTmp = new int[len];
	this.nextECsTmp = new int[len];
	
	for(int i=0; i < edgeList.size(); i++){
	    this.edgeIndicies[i] = i;
	    Edge e = edgeList.get(i);
	    curIndicies[i] = e.nthInterval(0).getThreadIndex();
	    nextECs[i] = 1; //nextEC is alwasy initialized to 1 since the firstEC starts at 0.
	}
    }
    
    public String printCurIndicies(){
	StringBuffer bf = new StringBuffer();
	for(int i=0; i<curIndicies.length; i++){
	    bf.append("\t" + curIndicies[i]);
	}
	return bf.toString();
    }

    public int[] getEdgeIndicies(){
	return this.edgeIndicies;
    }

    public int[] getCurIndicies(){
	return this.curIndicies;
    }

    public Edge getEdgeForIndex(ArrayList<Edge> edgeList, int i){
	return edgeList.get(edgeIndicies[i]);
    }

    public int getECForEdgeForIndex(int i){
	return (this.nextECs[i] - 1);
    }

    public int length(){
	return this.len;
    }
    
    public boolean feedNext(ArrayList<Edge> edgeList){
	//System.err.println("FEEDNEXT");
	return this.feedNext(0, edgeList);
    }
    //returns false if no more to feed
    private boolean feedNext(int i, ArrayList<Edge> el){
	if( this.nextECs[i] >= el.get(edgeIndicies[i]).getMultiplicity() )
	    return false;
	else{
	    curIndicies[i] = el.get(edgeIndicies[i]).nthInterval(this.nextECs[i]).getThreadIndex();
	    this.nextECs[i]++;
	    return true;
	}
    }

    /* feed all after sorted */
    public boolean feedAll(ArrayList<Edge> el){
	for(int i=0; i<len;i++){
	    if(!feedNext(i,el))
		return false;
	}
	return true;
    }

    /* this sorts 3 set of indicies based on curIndicies, which is threadIndex for current intervals for each edge*/
    public void sort(){
	this.curIndiciesTmp = new int[len];
	this.edgeIndiciesTmp = new int[len];
	this.nextECsTmp = new int[len];
	
	int lowerIndex = 0;
	int higherIndex = len - 1;
	
	this.mergesort(lowerIndex, higherIndex);
    }

    private void mergesort(int low, int high){
	if(low < high){
	    int middle = low + (high - low) / 2;
	    mergesort(low, middle);
	    mergesort(middle + 1, high);
	    merge(low, middle, high);
	}
    }

    private void copy2Tmp(int i){
	this.curIndiciesTmp[i] = this.curIndicies[i];
	this.edgeIndiciesTmp[i] = this.edgeIndicies[i];
	this.nextECsTmp[i] = this.nextECs[i];
    }

    private void setFromTmp( int origIndex, int tmpIndex){
	this.curIndicies[origIndex] = this.curIndiciesTmp[tmpIndex];
	this.edgeIndicies[origIndex] = this.edgeIndiciesTmp[tmpIndex];
	this.nextECs[origIndex] = this.nextECsTmp[tmpIndex];
    }

    private void merge(int low, int middle, int high){
	for (int i = low; i <= high; i++) {
	    this.copy2Tmp(i);
	}
	int i = low;
	int j = middle + 1;
	int k = low;
	while (i <= middle && j <= high) {
	    if (curIndiciesTmp[i] <= curIndiciesTmp[j]) {
		this.setFromTmp( k , i );
		i++;
		} else {
		this.setFromTmp( k , j );
		j++;
	    }
	    k++;
	}
	while (i <= middle) {
	    this.setFromTmp( k , i );
	    k++;
	    i++;
	}
    }

}
