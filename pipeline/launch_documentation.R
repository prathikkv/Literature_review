#!/usr/bin/env Rscript
#' Documentation Viewer Launcher
#' 
#' Simple launcher to open interactive documentation in default browser
#' 
#' @author Bioinformatics Pipeline Development Team
#' @version 1.0.0

#' Launch Interactive Documentation
#'
#' Opens the interactive HTML documentation in default web browser
#' @param doc_file Path to HTML documentation file
#' @return TRUE if successful
launch_documentation <- function(doc_file = "Interactive_Technical_Documentation.html") {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║               DOCUMENTATION VIEWER LAUNCHER                  ║\n")
  cat("║               Opening Interactive Documentation              ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Check if documentation file exists
  if (!file.exists(doc_file)) {
    cat("❌ ERROR: Documentation file not found:", doc_file, "\n")
    cat("💡 Run 'Rscript generate_interactive_documentation.R' first\n\n")
    return(FALSE)
  }
  
  cat("📄 Documentation file:", doc_file, "\n")
  
  # Get file info
  file_info <- file.info(doc_file)
  file_size_kb <- round(file_info$size / 1024, 1)
  last_modified <- format(file_info$mtime, "%Y-%m-%d %H:%M:%S")
  
  cat("📏 File size:", file_size_kb, "KB\n")
  cat("🕒 Last modified:", last_modified, "\n\n")
  
  # Get full path for browser
  full_path <- normalizePath(doc_file, winslash = "/")
  browser_url <- paste0("file://", full_path)
  
  cat("🌐 Opening in browser:", browser_url, "\n")
  
  # Open in default browser (cross-platform)
  tryCatch({
    if (Sys.info()["sysname"] == "Darwin") {
      # macOS
      system2("open", browser_url, wait = FALSE)
      cat("🍎 macOS: Opened with 'open' command\n")
    } else if (Sys.info()["sysname"] == "Windows") {
      # Windows
      system2("start", browser_url, wait = FALSE)
      cat("🪟 Windows: Opened with 'start' command\n")
    } else {
      # Linux and others
      system2("xdg-open", browser_url, wait = FALSE)
      cat("🐧 Linux: Opened with 'xdg-open' command\n")
    }
    
    cat("\n🎉 SUCCESS: Documentation opened in browser!\n")
    cat("\n📋 Interactive Features Available:\n")
    cat("  • 📊 70+ Interactive Mermaid Flowcharts\n")
    cat("  • 📋 Table of Contents Sidebar (Ctrl/Cmd+T)\n")
    cat("  • 📥 Download individual flowcharts as SVG\n")
    cat("  • 🔍 Zoom functionality for detailed viewing\n")
    cat("  • 📱 Responsive design for all screen sizes\n")
    cat("  • 🖨️  Print-friendly formatting\n")
    cat("  • ⌨️  Keyboard shortcuts (Escape to close TOC)\n")
    
    cat("\n💡 Pro Tips:\n")
    cat("  • Use the 'Table of Contents' button to navigate quickly\n")
    cat("  • Download flowcharts as SVG for presentations\n")
    cat("  • Zoom flowcharts for detailed analysis\n")
    cat("  • Print or save as PDF for offline viewing\n\n")
    
    return(TRUE)
    
  }, error = function(e) {
    cat("❌ ERROR opening browser:", e$message, "\n")
    cat("💡 Manually open this URL in your browser:\n")
    cat("   ", browser_url, "\n\n")
    return(FALSE)
  })
}

#' Quick Documentation Update and Launch
#'
#' Regenerates documentation and opens in browser
#' @return TRUE if successful
update_and_launch <- function() {
  cat("🔄 Updating documentation...\n")
  
  # Source and run the generator
  if (file.exists("generate_interactive_documentation.R")) {
    source("generate_interactive_documentation.R")
    result <- generate_interactive_documentation()
    
    if (result) {
      cat("✅ Documentation updated successfully\n")
      return(launch_documentation())
    } else {
      cat("❌ Documentation update failed\n")
      return(FALSE)
    }
  } else {
    cat("❌ ERROR: generate_interactive_documentation.R not found\n")
    return(FALSE)
  }
}

#' Show Documentation Status
#'
#' Shows current status of documentation files
show_status <- function() {
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║                  DOCUMENTATION STATUS                        ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Check source markdown
  md_file <- "Technical_Documentation_CAMK2D_Pipeline.md"
  html_file <- "Interactive_Technical_Documentation.html"
  generator_file <- "generate_interactive_documentation.R"
  
  files_to_check <- list(
    "📄 Source Markdown" = md_file,
    "🌐 Interactive HTML" = html_file,
    "🔧 Generator Script" = generator_file
  )
  
  for (name in names(files_to_check)) {
    file_path <- files_to_check[[name]]
    
    if (file.exists(file_path)) {
      file_info <- file.info(file_path)
      size_kb <- round(file_info$size / 1024, 1)
      modified <- format(file_info$mtime, "%Y-%m-%d %H:%M")
      cat("✅", name, ":", file_path, paste0("(", size_kb, " KB, ", modified, ")"), "\n")
    } else {
      cat("❌", name, ": NOT FOUND -", file_path, "\n")
    }
  }
  
  # Count flowcharts if markdown exists
  if (file.exists(md_file)) {
    content <- readLines(md_file, warn = FALSE)
    mermaid_count <- sum(grepl("```mermaid", content))
    cat("\n📊 Flowcharts available:", mermaid_count, "\n")
  }
  
  cat("\n")
}

# Main execution logic
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    # Default: just launch
    launch_documentation()
  } else if (args[1] == "update") {
    # Update and launch
    update_and_launch()
  } else if (args[1] == "status") {
    # Show status
    show_status()
  } else {
    cat("Usage:\n")
    cat("  Rscript launch_documentation.R         # Launch documentation\n")
    cat("  Rscript launch_documentation.R update  # Update and launch\n")
    cat("  Rscript launch_documentation.R status  # Show status\n\n")
  }
} else {
  # If running interactively, show status
  show_status()
}

cat("✅ Documentation Viewer Launcher loaded successfully\n\n")