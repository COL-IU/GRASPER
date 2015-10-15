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

class DRWFPosComp implements Comparator<DRWFCompareIndex>{
    public int compare(DRWFCompareIndex i1, DRWFCompareIndex i2){
	if(i1.getPosition() < i2.getPosition())
	    return -1;
	else if(i1.getPosition() == i2.getPosition())
	    return 0;
	else
	    return 1;
    }
}

//DirectionalRangeWithinFacingCompareIndex
public class DRWFCompareIndex{// implements Comparable<DRWFCompareIndex>{
    private int tDR;
    private int oDR;
    private boolean tDRFirst;
    private int size;
    private int position;
    //0: T-O
    //1: T-OP
    //2: TP-O
    //3: TP-OP
    private int pt;//pairingType 
    
    //if is WithinFacing is false, it's follower.
    public DRWFCompareIndex(DirectionalRange t, int ti, DirectionalRange o, int oi, boolean tFirst, int pairingType){
	this.tDR = ti;
	this.oDR = oi;
	this.tDRFirst = tFirst;
	// ---this---> <----other---
	if(tFirst){
	    this.size = t.getMax() - o.getMin() + 1;
	    this.position = t.getMax() + (o.getMin() - t.getMax())/2;
	}
	// --other---> <----this----
	else{
	    this.size = t.getMax() - o.getMin() + 1;
	    this.position = o.getMax() + (t.getMin() - o.getMax())/2;
	}
	this.pt = pairingType;
    }
    
    public int pairingType(){
	return this.pt;
    }

    public int tDRIndex(){
	return this.tDR;
    }
    public int oDRIndex(){
	return this.oDR;
    }
    public boolean tDRFirst(){
	return this.tDRFirst;
    }
    public int getSize(){
	return this.size;
    }

    public int getPosition(){
	return this.position;
    }

    
    //returns a list of directional range
    public DirectionalRange[] isInvertionPair(DRWFCompareIndex other, HalfCluster t, HalfCluster o){

	/*boolean debug = false;
	if( ( (t.getId() == 4221 || t.getPair().getId() == 4222) && (o.getId() == 4191 || o.getPair().getId() == 4192) )
	    ||
	    ( (o.getId() == 4221 || o.getPair().getId() == 4222) && (t.getId() == 4191 || t.getPair().getId() == 4192) )
	    )
	    debug = true;
	*/
	int tp = this.pt;
	int op = other.pt;
	int tPos = this.position;
	int oPos = other.getPosition();
	int iSize = oPos - tPos;
	boolean tFirst = true;
	if(iSize < 0){
	    iSize = 0 - iSize;
	    tFirst = false;
	}

	//if(debug)
	//  System.err.println("[4221] iSize :\t" + iSize);
	
	//first four is t,o,tp,op in correct order, last cell is for type 1 through 8
	DirectionalRange[] drs = new DirectionalRange[5];

	if(iSize >= Constants.MIN_INVERSION_SIZE 
	   && iSize <= Constants.MAX_INVERSION_SIZE){
	    
	    if(tp == 0 && op == 3){
		//     --T-><-O--    --TP-><-OP--    (I-1)
		// OR  --TP-><-OP--   --T-><-O---    (I-2)
		if(this.tDRFirst && other.tDRFirst()){
	
		    if(tFirst){
			drs[0] = t.getNthDirectionalRange(this.tDRIndex());
			drs[1] = o.getNthDirectionalRange(this.oDRIndex());
			drs[2] = t.getPair().getNthDirectionalRange(other.tDRIndex());
			drs[3] = o.getPair().getNthDirectionalRange(other.oDRIndex());
			drs[4] = new DirectionalRange(1,1);
		    }else{
			drs[0] = t.getPair().getNthDirectionalRange(other.tDRIndex());
			drs[1] = o.getPair().getNthDirectionalRange(other.oDRIndex());
			drs[2] = t.getNthDirectionalRange(this.tDRIndex());
			drs[3] = o.getNthDirectionalRange(this.oDRIndex());
			drs[4] = new DirectionalRange(2,2);
		    }		
		}
		//     --O-><-T--    --OP-><-TP--   (I-3)
		// OR  --OP-><-TP--   --O-><-T---   (I-4)
		else if(!this.tDRFirst && !other.tDRFirst()){
		    if(tFirst){
			drs[0] = o.getNthDirectionalRange(this.oDRIndex());
			drs[1] = t.getNthDirectionalRange(this.tDRIndex());
			drs[2] = o.getPair().getNthDirectionalRange(other.oDRIndex());
			drs[3] = t.getPair().getNthDirectionalRange(other.tDRIndex());
			drs[4] = new DirectionalRange(3,3);
		    }else{
			drs[0] = o.getPair().getNthDirectionalRange(other.oDRIndex());
			drs[1] = t.getPair().getNthDirectionalRange(other.tDRIndex());
			drs[2] = o.getNthDirectionalRange(this.oDRIndex());
			drs[3] = t.getNthDirectionalRange(this.tDRIndex());
			drs[4] = new DirectionalRange(4,4);
		    }		
		}else{
		    return null;
		}
	    }
	    else if(tp == 1 && op == 2){
		//     --T-><-OP--    --TP-><-O--    (I-5)
		// OR  --TP-><-O--   --T-><-OP---    (I-6)
		if(this.tDRFirst && other.tDRFirst()){
		    if(tFirst){
			drs[0] = t.getNthDirectionalRange(this.tDRIndex());
			drs[1] = o.getPair().getNthDirectionalRange(this.oDRIndex());
			drs[2] = t.getPair().getNthDirectionalRange(other.tDRIndex());
			drs[3] = o.getNthDirectionalRange(other.oDRIndex());
			drs[4] = new DirectionalRange(5,5);
		    }else{
			drs[0] = t.getPair().getNthDirectionalRange(other.tDRIndex());
			drs[1] = o.getNthDirectionalRange(other.oDRIndex());
			drs[2] = t.getNthDirectionalRange(this.tDRIndex());
			drs[3] = o.getPair().getNthDirectionalRange(this.oDRIndex());
			drs[4] = new DirectionalRange(6,6);
		    }
		}
		//     --OP-><-T--    --O-><-TP--    (I-7)
		// OR  --O-><-TP--   --OP-><-T---    (I-8)
		else if(!this.tDRFirst && !other.tDRFirst()){
		    if(tFirst){
			drs[0] = o.getPair().getNthDirectionalRange(this.oDRIndex());
			drs[1] = t.getNthDirectionalRange(this.tDRIndex());
			drs[2] = o.getNthDirectionalRange(other.oDRIndex());
			drs[3] = t.getPair().getNthDirectionalRange(other.tDRIndex());
			drs[4] = new DirectionalRange(7,7);
		    }else{
			drs[0] = o.getNthDirectionalRange(other.oDRIndex());
			drs[1] = t.getPair().getNthDirectionalRange(other.tDRIndex());
			drs[2] = o.getPair().getNthDirectionalRange(this.oDRIndex());
			drs[3] = t.getNthDirectionalRange(this.tDRIndex());
			drs[4] = new DirectionalRange(8,8);
		    }
		}else{
		    return null;
		}
	    }
	    /*if(debug){
		for(int i=0;i<4;i++)
		    System.err.println("DEBUG[4221]:drs["+i+"]:\t"+drs[i].toString());
		    }*/
		

	    int padding = drs[0].getPaddingWhenTestingisWithinFacingFwd(drs[1], drs[2], drs[3], true);
	    /*if(debug){
		System.err.println("DEBUG[4221]:\t" + padding);
		}*/
	    if(drs[0].isWithinFacingFwdINVERSION(drs[1], padding) && drs[2].isWithinFacingFwdINVERSION(drs[3], padding)){
		//if(debug)
		//   System.err.println("DEBUG[4221]:\tSUCCESS" );
		return drs;
	    }else
		return null;
	}
	return null;
    }
}
