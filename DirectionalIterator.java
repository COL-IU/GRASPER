import java.util.*;

//
//auxilary class to combine fwd and reverse iteration into single operation.
//
public class DirectionalIterator<E>{

    private ListIterator<E> itr;
    private boolean fwd;
    private boolean throwIllegalStateException;

    public DirectionalIterator(LinkedList<E> list, boolean fwd){
	this.fwd = fwd;

	if(fwd)
	    this.itr = list.listIterator();
	else // if reverse, then we need to traverse from the end of the list.
	    this.itr = list.listIterator(list.size());
    }
    
    public boolean hasNext(){
	if(fwd)
	    return this.itr.hasNext();
	else
	    return this.itr.hasPrevious();
    }
    
    public void add(E e){
	if(fwd)
	    this.itr.add(e);
	else{
	    this.itr.add(e);
	    this.itr.previous(); // need to put the cursor past the newly inserted element e.
	    this.throwIllegalStateException = true;
	}
    }

    public void addWhenNotNull(E e){
	this.previous();
	this.add(e);
	this.next();
	this.throwIllegalStateException = true;
    }
    
    public E next(){
	this.throwIllegalStateException = false;
	if(fwd)
	    return this.itr.next();
	else
	    return this.itr.previous();
	
    }
    
    public E previous(){
	this.throwIllegalStateException = false;
	if(fwd)
	    return this.itr.previous();
	else
	    return this.itr.next();
    }
    
    public void remove(){
	if(throwIllegalStateException)
	    throw new IllegalStateException("remove() called after remove or add have been called");
	else{
	    this.itr.remove();
	    this.throwIllegalStateException = true;
	}

    }
    
    public void set(E e){
	if(throwIllegalStateException)
	    throw new IllegalStateException("set(E e) called after remove or add have been called");
	else
	    this.itr.set(e);
    }
}
