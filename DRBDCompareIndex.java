import java.util.*;

class SizeComp implements Comparator<DRBDCompareIndex>{
    public int compare(DRBDCompareIndex i1, DRBDCompareIndex i2){
	if(i1.getSize() < i2.getSize())
	    return -1;
	else if(i1.getSize() == i2.getSize())
	    return 0;
	else
	    return 1;
    }
}

class PositionComp implements Comparator<DRBDCompareIndex>{
    public int compare(DRBDCompareIndex i1, DRBDCompareIndex i2){
	if(i1.getPosition() < i2.getPosition())
	    return -1;
	else if(i1.getPosition() == i2.getPosition())
	    return 0;
	else
	    return 1;
    }
}

//DirectionalRangeBoundaryDefinerCompareIndex
public class DRBDCompareIndex{// implements Comparable<DRBDCompareIndex>{
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
    
    //if is BoundaryDefiner is false, it's follower.
    public DRBDCompareIndex(DirectionalRange t, int ti, DirectionalRange o, int oi, boolean tFirst, int pairingType){
	this.tDR = ti;
	this.oDR = oi;
	this.tDRFirst = tFirst;
	//<---this---   ----other--->
	if(tFirst){
	    this.size = o.getMax() - t.getMin() + 1;
	    this.position = t.getMin();
	}
	//<--other---   ----this---->
	else{
	    this.size = t.getMax() - o.getMin() + 1;
	    this.position = o.getMin();
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

    public int getInvertedPT(){
	if(this.pt == 0)
	    return 3;
	else if(this.pt == 3)
	    return 0;
	else if(this.pt == 1)
	    return 2;
	else if(this.pt == 2)
	    return 1;
	System.err.println("INVALID pairingType: " + this.pt);
	System.exit(0);
	return -1;
    }

    /*
    public int compareTo(DRBDCompareIndex other){
	//diff = 0 : sameSize
	//diff > 0 : this is bigger;
	//diff < 0 : other is bigger;
	int diff = this.size - other.getSize();
	if(diff > Constants.MAX_D)
	    return 1;
	else if( diff > (0-Constants.MAX_D)){
	    int posDiff = this.position - other.getPosition();
	    if(posDiff > 0)
		return 1;
	    else if(posDiff == 0)
		return 0;
	    else
		return -1;
	}
	    
	}*/
}
