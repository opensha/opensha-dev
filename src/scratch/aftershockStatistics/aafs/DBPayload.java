package scratch.aftershockStatistics.aafs;

import scratch.aftershockStatistics.util.MarshalReader;
import scratch.aftershockStatistics.util.MarshalWriter;
import scratch.aftershockStatistics.util.MarshalException;

import scratch.aftershockStatistics.aafs.entity.PendingTask;
import scratch.aftershockStatistics.aafs.entity.LogEntry;
import scratch.aftershockStatistics.aafs.entity.CatalogSnapshot;
import scratch.aftershockStatistics.aafs.entity.TimelineEntry;


/**
 * Common base class for all payloads in the database.
 * Author: Michael Barall 05/05/2018.
 *
 * The design assumption is that the type of object can be determined from
 * the database record contents (opcode or actcode), and hence polymorphic
 * marshal/unmarshal is not required.
 */
public class DBPayload {

	//----- Constants and variables -----




	//----- Construction -----

	/**
	 * Default constructor does nothing.
	 */
	public DBPayload () {}




	//----- Marshaling -----

	// Marshal version number.

	private static final int MARSHAL_VER_1 = 25001;

	private static final String M_VERSION_NAME = "DBPayload";

	// Marshal object, internal.

	protected void do_marshal (MarshalWriter writer) {

		// Version

		writer.marshalInt (M_VERSION_NAME, MARSHAL_VER_1);

		// Contents

	
		return;
	}

	// Unmarshal object, internal.

	protected void do_umarshal (MarshalReader reader) {
	
		// Version

		int ver = reader.unmarshalInt (M_VERSION_NAME, MARSHAL_VER_1, MARSHAL_VER_1);

		// Contents


		return;
	}

	// Marshal object.

	public void marshal (MarshalWriter writer, String name) {
		writer.marshalMapBegin (name);
		do_marshal (writer);
		writer.marshalMapEnd ();
		return;
	}

	// Unmarshal object.

	public DBPayload unmarshal (MarshalReader reader, String name) {
		reader.unmarshalMapBegin (name);
		do_umarshal (reader);
		reader.unmarshalMapEnd ();
		return this;
	}

	// Marshal object, for a pending task.

	public MarshalWriter marshal_task () {
		MarshalWriter writer = PendingTask.begin_details();
		marshal (writer, null);
		return writer;
	}

	// Unmarshal object, for a pending task.

	public DBPayload unmarshal_task (PendingTask ptask) {
		try {
			unmarshal (ptask.get_details(), null);
		} catch (Exception e) {
			throw new DBCorruptException("Error unmarshaling pending task payload\n" + ptask.toString() + "\nDump:\n" + ptask.dump_details(), e);
		}
		return this;
	}

	// Marshal object, for a log entry.

	public MarshalWriter marshal_log () {
		MarshalWriter writer = LogEntry.begin_details();
		marshal (writer, null);
		return writer;
	}

	// Unmarshal object, for a log entry.

	public DBPayload unmarshal_log (LogEntry lentry) {
		try {
			unmarshal (lentry.get_details(), null);
		} catch (Exception e) {
			throw new DBCorruptException("Error unmarshaling log entry payload\n" + lentry.toString() + "\nDump:\n" + lentry.dump_details(), e);
		}
		return this;
	}

	// Marshal object, for a timeline entry.

	public MarshalWriter marshal_timeline () {
		MarshalWriter writer = TimelineEntry.begin_details();
		marshal (writer, null);
		return writer;
	}

	// Unmarshal object, for a timeline entry.

	public DBPayload unmarshal_timeline (TimelineEntry tentry) {
		try {
			unmarshal (tentry.get_details(), null);
		} catch (Exception e) {
			throw new DBCorruptException("Error unmarshaling timeline entry payload\n" + tentry.toString() + "\nDump:\n" + tentry.dump_details(), e);
		}
		return this;
	}

}
