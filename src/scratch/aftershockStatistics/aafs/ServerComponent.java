package scratch.aftershockStatistics.aafs;

/**
 * Common base class for server components.
 * Author: Michael Barall 06/23/2018.
 */
public class ServerComponent {


	//----- Constant definitions -----
	//
	// These constants are known to every server component.




	// Opcodes.

	public static final int OPCODE_MIN = 1;					// Minimum allowed opcode
	public static final int OPCODE_NO_OP = 1;				// No operation
	public static final int OPCODE_SHUTDOWN = 2;			// Shut down the server
	public static final int OPCODE_CON_MESSAGE = 3;			// Write message to console
	public static final int OPCODE_GEN_FORECAST = 4;		// Generate forecast
	public static final int OPCODE_GEN_PDL_REPORT = 5;		// Generate delayed forecast report to PDL
	public static final int OPCODE_GEN_EXPIRE = 6;			// Stop generation, mark the timeline expired
	public static final int OPCODE_INTAKE_SYNC = 7;			// Intake an event, from sync
	public static final int OPCODE_INTAKE_PDL = 8;			// Intake an event, from PDL
	public static final int OPCODE_ANALYST_INTERVENE = 9;	// Analyst intervention
	public static final int OPCODE_UNKNOWN = 10;			// Unknown operation
	public static final int OPCODE_ALIAS_SPLIT = 11;		// Notification of alias timeline split
	public static final int OPCODE_ALIAS_STOP = 12;			// Notification of alias timeline stop
	public static final int OPCODE_ALIAS_REVIVE = 13;		// Notification of alias timeline revive
	public static final int OPCODE_INTAKE_POLL = 14;		// Intake an event, from poll
	public static final int OPCODE_POLL_COMCAT_RUN = 15;	// Run a poll of Comcat to discover events
	public static final int OPCODE_POLL_COMCAT_START = 16;	// Start polling Comcat to discover events
	public static final int OPCODE_POLL_COMCAT_STOP = 17;	// Stop polling Comcat to discover events
	public static final int OPCODE_MAX = 17;				// Maximum allowed opcode

	// Return a string describing an opcode.

	public String get_opcode_as_string (int x) {
		switch (x) {
		case OPCODE_NO_OP: return "OPCODE_NO_OP";
		case OPCODE_SHUTDOWN: return "OPCODE_SHUTDOWN";
		case OPCODE_CON_MESSAGE: return "OPCODE_CON_MESSAGE";
		case OPCODE_GEN_FORECAST: return "OPCODE_GEN_FORECAST";
		case OPCODE_GEN_PDL_REPORT: return "OPCODE_GEN_PDL_REPORT";
		case OPCODE_GEN_EXPIRE: return "OPCODE_GEN_EXPIRE";
		case OPCODE_INTAKE_SYNC: return "OPCODE_INTAKE_SYNC";
		case OPCODE_INTAKE_PDL: return "OPCODE_INTAKE_PDL";
		case OPCODE_ANALYST_INTERVENE: return "OPCODE_ANALYST_INTERVENE";
		case OPCODE_UNKNOWN: return "OPCODE_UNKNOWN";
		case OPCODE_ALIAS_SPLIT: return "OPCODE_ALIAS_SPLIT";
		case OPCODE_ALIAS_STOP: return "OPCODE_ALIAS_STOP";
		case OPCODE_ALIAS_REVIVE: return "OPCODE_ALIAS_REVIVE";
		case OPCODE_INTAKE_POLL: return "OPCODE_INTAKE_POLL";
		case OPCODE_POLL_COMCAT_RUN: return "OPCODE_POLL_COMCAT_RUN";
		case OPCODE_POLL_COMCAT_START: return "OPCODE_POLL_COMCAT_START";
		case OPCODE_POLL_COMCAT_STOP: return "OPCODE_POLL_COMCAT_STOP";
		}
		return "OPCODE_INVALID(" + x + ")";
	}




	// Special execution times.

	public static final long EXEC_TIME_ACTIVE = 0L;						// Task is active
	public static final long EXEC_TIME_MIN_WAITING = 1L;				// Minimum execution time for waiting tasks
	public static final long EXEC_TIME_MIN_PROMPT = 10000000L;			// Minimum execution time for prompt tasks, 10^7 ~ 2.8 hours
	public static final long EXEC_TIME_MAX_PROMPT = 19000000L;			// Maximum execution time for prompt tasks
	public static final long EXEC_TIME_SHUTDOWN = 20000000L;			// Execution time for shutdown task, 2*10^7 ~ 5.6 hours
	public static final long EXEC_TIME_FAR_FUTURE = 1000000000000000L;	// 10^15 ~ 30,000 years




	// Result codes.

	public static final int RESCODE_MIN_NORMAL = 1;					// Minimum known normal result code
	public static final int RESCODE_SUCCESS = 1;					// Task completed successfully
	public static final int RESCODE_TASK_CORRUPT = 2;				// Task entry or payload was corrupted, task discarded
	public static final int RESCODE_TIMELINE_CORRUPT = 3;			// Timeline entry or payload was corrupted, task discarded
	public static final int RESCODE_TIMELINE_NOT_FOUND = 4;			// Timeline entry not found, task discarded
	public static final int RESCODE_TIMELINE_NOT_ACTIVE = 5;		// Timeline entry not active, task discarded
	public static final int RESCODE_TIMELINE_TASK_MISMATCH = 6;		// Timeline entry has lag values that do not match the forecast task
	public static final int RESCODE_TIMELINE_COMCAT_FAIL = 7;		// Timeline stopped due to ComCat failure
	public static final int RESCODE_TIMELINE_WITHDRAW = 8;			// Timeline stopped due to event not passing intake filter
	public static final int RESCODE_TIMELINE_FORESHOCK = 9;			// Timeline stopped because event was found to be a foreshock
	public static final int RESCODE_TIMELINE_NOT_PDL_PEND = 10;		// Timeline entry does not have a PDL report pending, task discarded
	public static final int RESCODE_TIMELINE_PDL_FAIL = 11;			// Timeline attempt to send PDL report failed, sending abandoned
	public static final int RESCODE_TIMELINE_EXISTS = 12;			// Timeline already exists, task discarded
	public static final int RESCODE_TASK_RETRY_SUCCESS = 13;		// Task completed on task dispatcher retry
	public static final int RESCODE_TIMELINE_STATE_UPDATE = 14;		// Timeline state was updated
	public static final int RESCODE_INTAKE_COMCAT_FAIL = 15;		// Event intake failed due to ComCat failure
	public static final int RESCODE_TIMELINE_ANALYST_SET = 16;		// Timeline analyst data was set
	public static final int RESCODE_TIMELINE_ANALYST_FAIL = 17;		// Timeline analyst intervention failed due to bad state
	public static final int RESCODE_TIMELINE_ANALYST_NONE = 18;		// Timeline analyst intervention not done
	public static final int RESCODE_ALIAS_TIMELINE_NOT_FOUND = 19;	// Timeline ID not found in the alias table
	public static final int RESCODE_ALIAS_STOPPED = 20;				// Timeline ID refers to a stopped timeline in the alias table
	public static final int RESCODE_ALIAS_NEW_EVENT = 21;			// Event ID does not appear in the alias table, can be a new timeline
	public static final int RESCODE_ALIAS_EVENT_NOT_IN_COMCAT = 22;	// Event ID is not known to Comcat, cannot query the alias table
	public static final int RESCODE_INTAKE_FILTERED = 23;			// Event intake dropped because event did not pass intake filter
	public static final int RESCODE_FORECAST_STALE = 24;			// Forecast skipped because it would be stale
	public static final int RESCODE_FORECAST_INTAKE_FILTERED = 25;	// Forecast skipped because event did not pass intake filter
	public static final int RESCODE_FORECAST_ANALYST_BLOCKED = 26;	// Forecast skipped because analyst blocked it
	public static final int RESCODE_FORECAST_SHADOWED = 27;			// Forecast skipped because the event is shadowed
	public static final int RESCODE_FORECAST_FORESHOCK = 28;		// Forecast skipped because event was found to be a foreshock
	public static final int RESCODE_INTAKE_AGED = 29;				// Event intake dropped because event was outside age range
	public static final int RESCODE_TIMELINE_EVENT_REMOVED = 30;	// Timeline stopped due to event removed or merged in Comcat
	public static final int RESCODE_MAX_NORMAL = 30;				// Maximum known normal result code

	public static final int RESCODE_DELETE = 101;					// Delete current task (without logging it)
	public static final int RESCODE_DELETE_TIMELINE_EXISTS = 102;	// Delete current task (without logging it), because timeline already exists
	public static final int RESCODE_DELETE_INTAKE_FILTERED = 103;	// Delete current task (without logging it), because event did not pass intake filter
	public static final int RESCODE_DELETE_INTAKE_AGED = 104;		// Delete current task (without logging it), because event intake was outside age range
	public static final int RESCODE_DELETE_TIMELINE_NO_ALIAS = 105;	// Delete current task (without logging it), because timeline has no active alias
	public static final int RESCODE_DELETE_TIMELINE_BAD_STATE = 106;	// Delete current task (without logging it), because timeline state does not permit operation
	public static final int RESCODE_DELETE_NOT_IN_COMCAT = 107;		// Delete current task (without logging it), because event is not in Comcat

	public static final int RESCODE_STAGE = 201;					// Stage current task (execute it again)
	public static final int RESCODE_STAGE_COMCAT_RETRY = 202;		// Stage current task (execute it again), to retry a failed Comcat operation
	public static final int RESCODE_STAGE_PDL_RETRY = 203;			// Stage current task (execute it again), to retry a failed PDL send
	public static final int RESCODE_STAGE_EVENT_ID = 204;			// Stage current task (execute it again), to insert event ID
	public static final int RESCODE_STAGE_TIMELINE_ID = 205;		// Stage current task (execute it again), to insert timeline ID
	public static final int RESCODE_STAGE_TOO_SOON = 206;			// Stage current task (execute it again), because it is too soon for the operation
	public static final int RESCODE_STAGE_REPEATING_TASK = 207;		// Stage current task (execute it again), because it is a repeating task


	// Return a string describing a result code.

	public String get_rescode_as_string (int x) {
		switch (x) {
		case RESCODE_SUCCESS: return "RESCODE_SUCCESS";
		case RESCODE_TASK_CORRUPT: return "RESCODE_TASK_CORRUPT";
		case RESCODE_TIMELINE_CORRUPT: return "RESCODE_TIMELINE_CORRUPT";
		case RESCODE_TIMELINE_NOT_FOUND: return "RESCODE_TIMELINE_NOT_FOUND";
		case RESCODE_TIMELINE_NOT_ACTIVE: return "RESCODE_TIMELINE_NOT_ACTIVE";
		case RESCODE_TIMELINE_TASK_MISMATCH: return "RESCODE_TIMELINE_TASK_MISMATCH";
		case RESCODE_TIMELINE_COMCAT_FAIL: return "RESCODE_TIMELINE_COMCAT_FAIL";
		case RESCODE_TIMELINE_WITHDRAW: return "RESCODE_TIMELINE_WITHDRAW";
		case RESCODE_TIMELINE_FORESHOCK: return "RESCODE_TIMELINE_FORESHOCK";
		case RESCODE_TIMELINE_NOT_PDL_PEND: return "RESCODE_TIMELINE_NOT_PDL_PEND";
		case RESCODE_TIMELINE_PDL_FAIL: return "RESCODE_TIMELINE_PDL_FAIL";
		case RESCODE_TIMELINE_EXISTS: return "RESCODE_TIMELINE_EXISTS";
		case RESCODE_TASK_RETRY_SUCCESS: return "RESCODE_TASK_RETRY_SUCCESS";
		case RESCODE_TIMELINE_STATE_UPDATE: return "RESCODE_TIMELINE_STATE_UPDATE";
		case RESCODE_INTAKE_COMCAT_FAIL: return "RESCODE_INTAKE_COMCAT_FAIL";
		case RESCODE_TIMELINE_ANALYST_SET: return "RESCODE_TIMELINE_ANALYST_SET";
		case RESCODE_TIMELINE_ANALYST_FAIL: return "RESCODE_TIMELINE_ANALYST_FAIL";
		case RESCODE_TIMELINE_ANALYST_NONE: return "RESCODE_TIMELINE_ANALYST_NONE";
		case RESCODE_ALIAS_TIMELINE_NOT_FOUND: return "RESCODE_ALIAS_TIMELINE_NOT_FOUND";
		case RESCODE_ALIAS_STOPPED: return "RESCODE_ALIAS_STOPPED";
		case RESCODE_ALIAS_NEW_EVENT: return "RESCODE_ALIAS_NEW_EVENT";
		case RESCODE_ALIAS_EVENT_NOT_IN_COMCAT: return "RESCODE_ALIAS_EVENT_NOT_IN_COMCAT";
		case RESCODE_INTAKE_FILTERED: return "RESCODE_INTAKE_FILTERED";
		case RESCODE_FORECAST_STALE: return "RESCODE_FORECAST_STALE";
		case RESCODE_FORECAST_INTAKE_FILTERED: return "RESCODE_FORECAST_INTAKE_FILTERED";
		case RESCODE_FORECAST_ANALYST_BLOCKED: return "RESCODE_FORECAST_ANALYST_BLOCKED";
		case RESCODE_FORECAST_SHADOWED: return "RESCODE_FORECAST_SHADOWED";
		case RESCODE_FORECAST_FORESHOCK: return "RESCODE_FORECAST_FORESHOCK";
		case RESCODE_INTAKE_AGED: return "RESCODE_INTAKE_AGED";
		case RESCODE_TIMELINE_EVENT_REMOVED: return "RESCODE_TIMELINE_EVENT_REMOVED";

		case RESCODE_DELETE: return "RESCODE_DELETE";
		case RESCODE_DELETE_TIMELINE_EXISTS: return "RESCODE_DELETE_TIMELINE_EXISTS";
		case RESCODE_DELETE_INTAKE_FILTERED: return "RESCODE_DELETE_INTAKE_FILTERED";
		case RESCODE_DELETE_INTAKE_AGED: return "RESCODE_DELETE_INTAKE_AGED";
		case RESCODE_DELETE_TIMELINE_NO_ALIAS: return "RESCODE_DELETE_TIMELINE_NO_ALIAS";
		case RESCODE_DELETE_TIMELINE_BAD_STATE: return "RESCODE_DELETE_TIMELINE_BAD_STATE";
		case RESCODE_DELETE_NOT_IN_COMCAT: return "RESCODE_DELETE_NOT_IN_COMCAT";

		case RESCODE_STAGE: return "RESCODE_STAGE";
		case RESCODE_STAGE_COMCAT_RETRY: return "RESCODE_STAGE_COMCAT_RETRY";
		case RESCODE_STAGE_PDL_RETRY: return "RESCODE_STAGE_PDL_RETRY";
		case RESCODE_STAGE_EVENT_ID: return "RESCODE_STAGE_EVENT_ID";
		case RESCODE_STAGE_TIMELINE_ID: return "RESCODE_STAGE_TIMELINE_ID";
		case RESCODE_STAGE_TOO_SOON: return "RESCODE_STAGE_TOO_SOON";
		case RESCODE_STAGE_REPEATING_TASK: return "RESCODE_STAGE_REPEATING_TASK";
		}
		return "RESCODE_INVALID(" + x + ")";
	}




	// Special event ids.

	public static final String EVID_SHUTDOWN = "===shutdown===";	// Shutdown task

	public static final String EVID_UNKNOWN = "===unknown===";		// Unknown event ID

	public static final String EVID_ERROR = "===error===";	// Error status

	public static final String EVID_ANALYST = "===analyst===";	// Used for analyst-selected shadowing

	public static final String EVID_POLL = "===poll===";	// Used for polling tasks




	// Special submit ids.

	public static final String SUBID_AAFS = "AAFS";		// Automatic system




	// Special durations, in milliseconds.

	public static final long DURATION_DAY = 86400000L;				// 1 day
	public static final long DURATION_YEAR = 31536000000L;			// 1 year (365 days)
	public static final long DURATION_HUGE = 1000000000000000L;		// 10^15 ~ 30,000 years




	// Logging.

	public static final String LOG_SEPARATOR_LINE = "------------------------------";	// Used to separate log entries




	//----- Component access -----


	// Server group.
	// This provides access to other server components.

	protected ServerGroup sg;




	//----- Construction -----


	// Default constructor.

	public ServerComponent () {
		this.sg = null;
	}


	// Set up this component by linking to the server group.
	// A subclass may override this to perform additional setup operations.

	public void setup (ServerGroup the_sg) {
		this.sg = the_sg;
		return;
	}

}
