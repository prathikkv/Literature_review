#!/usr/bin/env Rscript
#' Error Handling and Logging System
#' 
#' Comprehensive error handling, logging, and monitoring for bioinformatics pipeline
#' Provides structured error reporting, performance tracking, and audit trails
#' 
#' @author Claude Code Pipeline Standards
#' @version 1.0.0

# Load required libraries
suppressPackageStartupMessages({
  library(jsonlite)
  library(digest)
})

# Global error handler settings
.error_handler_settings <- list(
  log_level = "INFO",
  max_log_size_mb = 100,
  log_rotation_enabled = TRUE,
  console_output = TRUE,
  email_alerts_enabled = FALSE,
  error_report_path = "output/logs/error_reports"
)

#' Initialize Logging System
#'
#' Sets up logging directories and configuration
#' @param config Pipeline configuration
#' @param run_id Current run identifier
#' @return Logging configuration
initialize_logging_system <- function(config, run_id) {
  
  # Create log directories
  log_dirs <- list(
    base = "output/logs",
    current_run = file.path("output/runs", run_id, "logs"),
    error_reports = "output/logs/error_reports",
    performance = "output/logs/performance",
    security = "output/logs/security"
  )
  
  for (dir_path in log_dirs) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  # Initialize log files
  log_files <- list(
    main = file.path(log_dirs$current_run, "pipeline.log"),
    error = file.path(log_dirs$current_run, "errors.log"),
    performance = file.path(log_dirs$performance, paste0("performance_", run_id, ".log")),
    security = file.path(log_dirs$security, "security_audit.log"),
    step_timing = file.path(log_dirs$current_run, "step_timing.log")
  )
  
  # Write initial log entries
  for (log_file in log_files) {
    if (!file.exists(log_file)) {
      cat("# Pipeline Log Started:", as.character(Sys.time()), "\n", file = log_file)
    }
  }
  
  return(list(
    directories = log_dirs,
    files = log_files,
    initialized_at = Sys.time(),
    run_id = run_id
  ))
}

#' Log Message
#'
#' Central logging function with multiple severity levels
#' @param message Log message
#' @param level Log level ("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL")
#' @param category Category ("PIPELINE", "SECURITY", "PERFORMANCE", "DATA")
#' @param details Additional details (optional)
#' @param log_config Logging configuration
log_message <- function(message, 
                       level = "INFO", 
                       category = "PIPELINE",
                       details = NULL,
                       log_config = NULL) {
  
  # Create timestamp
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  # Format log entry
  log_entry <- paste(
    timestamp,
    level,
    category,
    message,
    sep = " | "
  )
  
  # Add details if provided
  if (!is.null(details)) {
    log_entry <- paste(log_entry, "\nDetails:", jsonlite::toJSON(details, auto_unbox = TRUE))
  }
  
  # Console output
  if (.error_handler_settings$console_output) {
    # Color coding for console
    color_code <- switch(level,
      "DEBUG" = "\033[36m",    # Cyan
      "INFO" = "\033[32m",     # Green
      "WARNING" = "\033[33m",  # Yellow
      "ERROR" = "\033[31m",    # Red
      "CRITICAL" = "\033[35m", # Magenta
      ""
    )
    reset_code <- "\033[0m"
    
    cat(paste0(color_code, "[", level, "] ", message, reset_code, "\n"))
  }
  
  # Write to log files
  if (!is.null(log_config)) {
    # Main log
    cat(log_entry, "\n", file = log_config$files$main, append = TRUE)
    
    # Error log for errors and warnings
    if (level %in% c("ERROR", "CRITICAL", "WARNING")) {
      cat(log_entry, "\n", file = log_config$files$error, append = TRUE)
    }
    
    # Security log for security-related events
    if (category == "SECURITY") {
      cat(log_entry, "\n", file = log_config$files$security, append = TRUE)
    }
    
    # Performance log for performance events
    if (category == "PERFORMANCE") {
      cat(log_entry, "\n", file = log_config$files$performance, append = TRUE)
    }
  }
  
  # Check for log rotation
  if (!is.null(log_config)) {
    check_log_rotation(log_config)
  }
}

#' Enhanced Error Handler
#'
#' Comprehensive error handling with context capture
#' @param error Error object
#' @param context Additional context information
#' @param step_name Current pipeline step
#' @param recovery_action Suggested recovery action
#' @param log_config Logging configuration
#' @return Error report
handle_pipeline_error <- function(error, 
                                 context = NULL,
                                 step_name = "unknown",
                                 recovery_action = "Contact support",
                                 log_config = NULL) {
  
  # Capture system state
  system_state <- capture_system_state()
  
  # Create comprehensive error report
  error_report <- list(
    timestamp = as.character(Sys.time()),
    error_id = generate_error_id(),
    step_name = step_name,
    error_message = as.character(error$message),
    error_class = class(error),
    call_stack = capture_call_stack(),
    system_state = system_state,
    context = context,
    recovery_action = recovery_action,
    severity = assess_error_severity(error)
  )
  
  # Log the error
  log_message(
    message = paste("Pipeline Error in", step_name, ":", error$message),
    level = "ERROR",
    category = "PIPELINE",
    details = list(
      error_id = error_report$error_id,
      severity = error_report$severity,
      recovery = recovery_action
    ),
    log_config = log_config
  )
  
  # Save detailed error report
  save_error_report(error_report, log_config)
  
  # Check if error is critical
  if (error_report$severity == "CRITICAL") {
    handle_critical_error(error_report, log_config)
  }
  
  return(error_report)
}

#' Capture System State
#'
#' Captures current system state for debugging
#' @return System state information
capture_system_state <- function() {
  
  tryCatch({
    list(
      memory_usage = capture_memory_usage(),
      disk_space = capture_disk_space(),
      r_session = capture_r_session_info(),
      environment_variables = capture_env_vars(),
      working_directory = getwd(),
      loaded_packages = capture_loaded_packages(),
      system_load = capture_system_load()
    )
  }, error = function(e) {
    list(error = "Failed to capture system state", details = e$message)
  })
}

#' Capture Memory Usage
#'
#' Gets current memory usage statistics
#' @return Memory usage information
capture_memory_usage <- function() {
  
  tryCatch({
    if (.Platform$OS.type == "unix") {
      # Unix-like systems
      gc_info <- gc()
      list(
        r_memory_mb = sum(gc_info[, "used"]) * 0.048576, # Rough conversion
        gc_level0 = gc_info[1, "used"],
        gc_level1 = gc_info[2, "used"],
        system_memory = "unavailable_unix"
      )
    } else {
      # Windows systems
      list(
        r_memory_mb = sum(gc()[, "used"]) * 0.048576,
        system_memory = "unavailable_windows"
      )
    }
  }, error = function(e) {
    list(error = "Memory capture failed", details = e$message)
  })
}

#' Capture Disk Space
#'
#' Gets available disk space
#' @return Disk space information
capture_disk_space <- function() {
  
  tryCatch({
    if (.Platform$OS.type == "unix") {
      # Try to get disk space on Unix
      df_output <- system("df -h .", intern = TRUE)
      list(
        command_output = df_output,
        working_dir = getwd()
      )
    } else {
      list(
        status = "disk_space_check_unavailable_windows",
        working_dir = getwd()
      )
    }
  }, error = function(e) {
    list(error = "Disk space check failed", details = e$message)
  })
}

#' Capture R Session Info
#'
#' Gets R session information
#' @return R session details
capture_r_session_info <- function() {
  
  tryCatch({
    session <- sessionInfo()
    list(
      r_version = paste(session$R.version$major, session$R.version$minor, sep = "."),
      platform = session$platform,
      locale = session$locale,
      base_packages = names(session$basePkgs),
      attached_packages = if(!is.null(session$otherPkgs)) names(session$otherPkgs) else character(0)
    )
  }, error = function(e) {
    list(error = "R session capture failed", details = e$message)
  })
}

#' Capture Environment Variables
#'
#' Gets relevant environment variables (safely)
#' @return Environment variables
capture_env_vars <- function() {
  
  safe_vars <- c("R_HOME", "R_LIBS", "PATH", "HOME", "USER", "LANG")
  
  env_vars <- list()
  for (var in safe_vars) {
    env_vars[[var]] <- Sys.getenv(var, unset = "not_set")
  }
  
  return(env_vars)
}

#' Capture Loaded Packages
#'
#' Gets information about loaded packages
#' @return Package information
capture_loaded_packages <- function() {
  
  tryCatch({
    attached <- search()
    loaded <- loadedNamespaces()
    
    list(
      attached_packages = attached[grepl("package:", attached)],
      loaded_namespaces = loaded,
      package_count = length(loaded)
    )
  }, error = function(e) {
    list(error = "Package capture failed", details = e$message)
  })
}

#' Capture System Load
#'
#' Gets basic system load information
#' @return System load details
capture_system_load <- function() {
  
  tryCatch({
    list(
      timestamp = as.character(Sys.time()),
      processes = length(list.files("/proc", pattern = "^[0-9]+$", full.names = FALSE)),
      r_process_id = Sys.getpid()
    )
  }, error = function(e) {
    list(timestamp = as.character(Sys.time()), error = "System load capture failed")
  })
}

#' Capture Call Stack
#'
#' Captures the current call stack for debugging
#' @return Call stack information
capture_call_stack <- function() {
  
  tryCatch({
    calls <- sys.calls()
    call_info <- list()
    
    for (i in 1:min(length(calls), 10)) {  # Limit to last 10 calls
      call_info[[i]] <- list(
        call_number = i,
        function_name = deparse(calls[[i]][[1]]),
        call_text = deparse(calls[[i]])
      )
    }
    
    return(call_info)
  }, error = function(e) {
    return(list(error = "Call stack capture failed"))
  })
}

#' Generate Error ID
#'
#' Creates unique error identifier
#' @return Unique error ID
generate_error_id <- function() {
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  random_part <- sample(10000:99999, 1)
  
  return(paste0("ERR_", timestamp, "_", random_part))
}

#' Assess Error Severity
#'
#' Determines error severity level
#' @param error Error object
#' @return Severity level
assess_error_severity <- function(error) {
  
  error_msg <- tolower(as.character(error$message))
  
  # Critical errors - system/data integrity issues
  if (any(grepl("file not found|permission denied|disk space|memory|corruption", error_msg))) {
    return("CRITICAL")
  }
  
  # High severity - analysis failures
  if (any(grepl("failed to|cannot|invalid|missing", error_msg))) {
    return("HIGH")
  }
  
  # Medium severity - warnings that became errors
  if (any(grepl("warning|deprecated|suggest", error_msg))) {
    return("MEDIUM")
  }
  
  # Default to medium
  return("MEDIUM")
}

#' Save Error Report
#'
#' Saves detailed error report to file
#' @param error_report Error report object
#' @param log_config Logging configuration
save_error_report <- function(error_report, log_config) {
  
  if (is.null(log_config)) {
    return(FALSE)
  }
  
  tryCatch({
    # Create error report filename
    report_file <- file.path(
      .error_handler_settings$error_report_path,
      paste0(error_report$error_id, "_error_report.json")
    )
    
    # Ensure directory exists
    dir.create(dirname(report_file), recursive = TRUE, showWarnings = FALSE)
    
    # Save as JSON
    writeLines(
      jsonlite::toJSON(error_report, pretty = TRUE, auto_unbox = TRUE),
      report_file
    )
    
    log_message(
      message = paste("Error report saved:", basename(report_file)),
      level = "INFO",
      category = "PIPELINE",
      log_config = log_config
    )
    
    return(TRUE)
    
  }, error = function(e) {
    cat("Failed to save error report:", e$message, "\n")
    return(FALSE)
  })
}

#' Handle Critical Error
#'
#' Special handling for critical errors
#' @param error_report Error report
#' @param log_config Logging configuration
handle_critical_error <- function(error_report, log_config) {
  
  log_message(
    message = "CRITICAL ERROR DETECTED - Pipeline execution halted",
    level = "CRITICAL",
    category = "PIPELINE",
    details = list(
      error_id = error_report$error_id,
      step = error_report$step_name,
      recovery = error_report$recovery_action
    ),
    log_config = log_config
  )
  
  # Create emergency backup of current state
  create_emergency_backup(log_config)
  
  # Additional notifications could be added here
  # (email, Slack, etc. if configured)
}

#' Create Emergency Backup
#'
#' Creates backup of current pipeline state during critical errors
#' @param log_config Logging configuration
create_emergency_backup <- function(log_config) {
  
  tryCatch({
    backup_dir <- file.path("output/emergency_backups", format(Sys.time(), "%Y%m%d_%H%M%S"))
    dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Backup key directories
    backup_targets <- c("output/current", "config.yml", "output/logs")
    
    for (target in backup_targets) {
      if (file.exists(target) || dir.exists(target)) {
        if (dir.exists(target)) {
          file.copy(target, backup_dir, recursive = TRUE)
        } else {
          file.copy(target, backup_dir)
        }
      }
    }
    
    log_message(
      message = paste("Emergency backup created:", backup_dir),
      level = "INFO",
      category = "PIPELINE",
      log_config = log_config
    )
    
  }, error = function(e) {
    cat("Failed to create emergency backup:", e$message, "\n")
  })
}

#' Check Log Rotation
#'
#' Manages log file rotation to prevent excessive disk usage
#' @param log_config Logging configuration
check_log_rotation <- function(log_config) {
  
  if (!.error_handler_settings$log_rotation_enabled) {
    return(invisible(NULL))
  }
  
  max_size_bytes <- .error_handler_settings$max_log_size_mb * 1024 * 1024
  
  for (log_file in log_config$files) {
    if (file.exists(log_file)) {
      file_size <- file.info(log_file)$size
      
      if (file_size > max_size_bytes) {
        rotate_log_file(log_file)
      }
    }
  }
}

#' Rotate Log File
#'
#' Rotates a log file when it gets too large
#' @param log_file Path to log file
rotate_log_file <- function(log_file) {
  
  tryCatch({
    # Create rotated filename
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    rotated_file <- paste0(log_file, ".", timestamp, ".old")
    
    # Move current log to rotated name
    file.rename(log_file, rotated_file)
    
    # Create new log file
    cat("# Log rotated at:", as.character(Sys.time()), "\n", file = log_file)
    
    cat("Log file rotated:", basename(log_file), "\n")
    
  }, error = function(e) {
    cat("Failed to rotate log file:", e$message, "\n")
  })
}

#' Performance Monitor
#'
#' Monitors and logs performance metrics
#' @param step_name Current step name
#' @param start_time Step start time
#' @param end_time Step end time
#' @param memory_before Memory before step
#' @param memory_after Memory after step
#' @param log_config Logging configuration
log_performance_metrics <- function(step_name,
                                  start_time,
                                  end_time = Sys.time(),
                                  memory_before = NULL,
                                  memory_after = NULL,
                                  log_config = NULL) {
  
  # Calculate metrics
  duration_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  performance_data <- list(
    step_name = step_name,
    start_time = as.character(start_time),
    end_time = as.character(end_time),
    duration_seconds = duration_seconds,
    duration_formatted = format_duration(duration_seconds),
    memory_before = memory_before,
    memory_after = memory_after
  )
  
  # Log performance
  log_message(
    message = paste("Step completed:", step_name, "in", format_duration(duration_seconds)),
    level = "INFO",
    category = "PERFORMANCE",
    details = performance_data,
    log_config = log_config
  )
  
  # Write to step timing log
  if (!is.null(log_config)) {
    timing_entry <- paste(
      as.character(start_time),
      step_name,
      duration_seconds,
      sep = ","
    )
    cat(timing_entry, "\n", file = log_config$files$step_timing, append = TRUE)
  }
  
  return(performance_data)
}

#' Format Duration
#'
#' Formats duration in human-readable format
#' @param seconds Duration in seconds
#' @return Formatted duration string
format_duration <- function(seconds) {
  
  if (seconds < 60) {
    return(paste(round(seconds, 2), "seconds"))
  } else if (seconds < 3600) {
    minutes <- floor(seconds / 60)
    remaining_seconds <- seconds %% 60
    return(paste(minutes, "minutes", round(remaining_seconds, 1), "seconds"))
  } else {
    hours <- floor(seconds / 3600)
    remaining_minutes <- floor((seconds %% 3600) / 60)
    return(paste(hours, "hours", remaining_minutes, "minutes"))
  }
}

#' Cleanup Logs
#'
#' Cleans up old log files to manage disk space
#' @param max_age_days Maximum age of logs to keep
#' @param log_base_dir Base log directory
cleanup_old_logs <- function(max_age_days = 30, log_base_dir = "output/logs") {
  
  if (!dir.exists(log_base_dir)) {
    return(invisible(NULL))
  }
  
  current_time <- Sys.time()
  log_files <- list.files(log_base_dir, pattern = "\\.log$", recursive = TRUE, full.names = TRUE)
  
  for (log_file in log_files) {
    file_age <- difftime(current_time, file.info(log_file)$mtime, units = "days")
    
    if (as.numeric(file_age) > max_age_days) {
      tryCatch({
        file.remove(log_file)
        cat("Removed old log file:", basename(log_file), "\n")
      }, error = function(e) {
        cat("Failed to remove old log file:", basename(log_file), e$message, "\n")
      })
    }
  }
}

# Module loading confirmation
cat("✅ Error Handling and Logging System loaded successfully\n")
cat("   Functions: initialize_logging_system(), log_message(),\n")
cat("             handle_pipeline_error(), log_performance_metrics()\n")
cat("   Features: Structured logging, error reporting, performance monitoring\n")
cat("   Version: 1.0.0 (Production Ready)\n\n")