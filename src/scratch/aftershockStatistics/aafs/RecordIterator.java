package scratch.aftershockStatistics.aafs;

import java.util.Iterator;
import java.util.NoSuchElementException;

import org.mongodb.morphia.query.MorphiaIterator;

/**
 * Iterator for iterating through records in a table (collection).
 * Author: Michael Barall 04/09/2018.
 *
 * This is used to iterate through database entries, when it may not be possible
 * to bring all the entries into memory as a List<T>.
 *
 * It implements Iterator<T> so that it is possible to step through the database
 * entries with standard iterator operations hasNext() and next().  The remove()
 * operation is passed thru to Morphia but I'm not sure if it works.
 *
 * It implements Iterable<T> so that it can be used in the foreach statement.
 *
 * It implements AuotCloseable so that it can be created in a try-with-resources
 * statement, which is recommended to ensure the Morphia iterator is closed.
 *
 * Only code very close to the database engine should create these objects or
 * access their contents.  All other code should treat these objects as opaque.
 */
public class RecordIterator<T> implements Iterator<T>, Iterable<T>, AutoCloseable {

	//----- Opaque contents -----

	// The Morphia iterator.

    private MorphiaIterator<T,T> morphia_iterator;

	// Constructor saves the Morphia iterator.

	public RecordIterator (MorphiaIterator<T,T> morphia_iterator) {
		this.morphia_iterator = morphia_iterator;
	}

	// Get the Morphia iterator.

	public MorphiaIterator<T,T> get_morphia_iterator () {
		return morphia_iterator;
	}

	//----- Implementation of Iterator<T> -----

	// Returns true if the iteration has more elements.
	// (In other words, returns true if next() would return an element rather than throwing an exception.)

	@Override
	public boolean hasNext() {
		return morphia_iterator.hasNext();
	}

	// Returns the next element in the iteration.
	// Throws NoSuchElementException if the iteration has no more elements.

	@Override
	public T next() {
		return morphia_iterator.next();
	}

	// Removes from the underlying collection the last element returned by this iterator (optional operation).
	// This method can be called only once per call to next().
	// The behavior of an iterator is unspecified if the underlying collection is modified while the iteration
	// is in progress in any way other than by calling this method.
	// Throws UnsupportedOperationException if the remove operation is not supported by this iterator.
	// Throws IllegalStateException if the next method has not yet been called, or the remove method
	// has already been called after the last call to the next method

	@Override
	public void remove() {
		morphia_iterator.remove();
		return;
	}

	//----- Implementation of Iterable<T> -----

	// Returns an iterator over a set of elements of type T.

	@Override
	public Iterator<T> iterator() {
		return this;
	}

	//----- Implementation of AutoCloseable -----

	// Closes this resource, relinquishing any underlying resources.
	// This method is invoked automatically on objects managed by the try-with-resources statement. 

	@Override
	public void close() {
		morphia_iterator.close();
		return;
	}

}
