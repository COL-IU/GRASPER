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

public class ActiveClusterList{

    private LinkedList<HalfCluster> clusterList;
    private LinkedList<Integer> positionReadRecentECforActiveClusterList;
    //private HashMap<Integer> clusterIDHash;

    public ActiveClusterList(){
	this.clusterList = new LinkedList<HalfCluster>();
	this.positionReadRecentECforActiveClusterList = new LinkedList<Integer>();
	//this.clusterIDHash = new HashMap<Integer>();
    }

    public boolean add(HalfCluster hc, int posReadCurEC){
	if(!this.clusterList.add(hc))
	    return false;
	if(!this.positionReadRecentECforActiveClusterList.add(new Integer(posReadCurEC))){
	    this.clusterList.removeLast();
	    return false;
	}
	return true;
    }

    public int size(){
	return this.clusterList.size();
    }
    
    public ActiveClusterListIterator listIterator(){
	return this.listIterator(true);
    }

    public ActiveClusterListIterator listIterator(boolean fwd){
	return new ActiveClusterListIterator(this.clusterList, this.positionReadRecentECforActiveClusterList, true);
    }

    /*
    public DirectionalIterator<HalfCluster> listIterator(boolean fwd){
	//return this.clusterList.listIterator(st);
	return new DirectionalIterator(this.clusterList, fwd);
    }
    
    public DirectionalIterator<Integer> listIteratorEC(boolean fwd){
	//return this.positionReadRecentECforActiveClusterList.listIterator(st);
	return new DirectionalIterator(this.positionReadRecentECforActiveClusterList, fwd);
	}*/


    /*
      public ListIterator<HalfCluster> listIterator(int st){
	return this.clusterList.listIterator(st);
    }
    
    public ListIterator<Integer> listIteratorEC(int st){
	return this.positionReadRecentECforActiveClusterList.listIterator(st);
    }
    */
    

}

class ActiveClusterListIterator{

    private DirectionalIterator<HalfCluster> cluIterator;
    private DirectionalIterator<Integer> positionReadRecentECIterator;

    ActiveClusterListIterator(LinkedList<HalfCluster> cluList, LinkedList<Integer> posReadRecentECforActiveCluList, boolean fwd){
	this.cluIterator = new DirectionalIterator<HalfCluster>(cluList, fwd);
	this.positionReadRecentECIterator = new DirectionalIterator<Integer>(posReadRecentECforActiveCluList, fwd);
    }
    
    public HalfCluster next(){
	this.positionReadRecentECIterator.next();
	return this.cluIterator.next();
    }

    public ClusterRecentECPair nextPair(){
	return new ClusterRecentECPair(this.cluIterator.next(), this.positionReadRecentECIterator.next());
    }
    
    public boolean hasNext(){
	return this.cluIterator.hasNext();
    }
    
    public void remove(){
	this.positionReadRecentECIterator.remove();
	this.cluIterator.remove();
    }
    
}

class ClusterRecentECPair{
    
    private HalfCluster cl;
    private int ec;
    
    public ClusterRecentECPair(HalfCluster cl, int ec){
	this.cl = cl;
	this.ec = ec;
    }
    
    public ClusterRecentECPair(HalfCluster cl, Integer ec){
	this(cl, ec.intValue());
    }
    
    public HalfCluster getClu(){
	return this.cl;
    }
    
    public int getRecentEC(){
	return this.ec;
    }
}
