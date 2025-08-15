#!/usr/bin/env Rscript
#' Interactive Documentation Generator
#' 
#' Converts Technical_Documentation_CAMK2D_Pipeline.md to interactive HTML
#' with live Mermaid.js flowchart rendering
#' 
#' @author Bioinformatics Pipeline Development Team
#' @version 1.0.0

# Load required libraries
suppressPackageStartupMessages({
  library(rmarkdown)
  library(knitr)
  library(htmltools)
  library(stringr)
})

#' Generate Interactive Documentation
#'
#' Converts markdown documentation to interactive HTML with Mermaid rendering
#' @param input_file Input markdown file path
#' @param output_file Output HTML file path
#' @param title Document title
#' @return TRUE if successful, FALSE otherwise
generate_interactive_documentation <- function(
  input_file = "Technical_Documentation_CAMK2D_Pipeline.md",
  output_file = "Interactive_Technical_Documentation.html",
  title = "CAMK2D Pipeline - Interactive Technical Documentation",
  include_results_summary = TRUE
) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║         INTERACTIVE DOCUMENTATION GENERATOR                  ║\n")
  cat("║         Converting Markdown to Interactive HTML              ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  cat("📄 Input file:", input_file, "\n")
  cat("🌐 Output file:", output_file, "\n")
  cat("📋 Title:", title, "\n\n")
  
  # Check if input file exists
  if (!file.exists(input_file)) {
    cat("❌ ERROR: Input file not found:", input_file, "\n")
    return(FALSE)
  }
  
  cat("📖 Reading markdown content...\n")
  
  # Read markdown content
  tryCatch({
    markdown_content <- readLines(input_file, warn = FALSE)
    total_lines <- length(markdown_content)
    cat("✅ Read", total_lines, "lines from input file\n")
    
    # Count Mermaid flowcharts
    mermaid_count <- sum(grepl("```mermaid", markdown_content))
    cat("📊 Found", mermaid_count, "Mermaid flowcharts\n\n")
    
  }, error = function(e) {
    cat("❌ ERROR reading input file:", e$message, "\n")
    return(FALSE)
  })
  
  cat("🔧 Processing markdown content...\n")
  
  # Create enhanced markdown content with HTML integration
  enhanced_content <- process_markdown_for_html(markdown_content, title)
  
  # Add analysis results summary if requested
  if (include_results_summary) {
    cat("📊 Adding analysis results summary...\n")
    results_summary <- generate_results_summary()
    enhanced_content <- c(enhanced_content[1:10], results_summary, enhanced_content[11:length(enhanced_content)])
  }
  
  cat("✅ Content processing complete\n\n")
  
  cat("🌐 Generating interactive HTML...\n")
  
  # Create temporary enhanced markdown file
  temp_md <- tempfile(fileext = ".md")
  writeLines(enhanced_content, temp_md)
  
  # Generate HTML with custom template
  tryCatch({
    # Create custom HTML output
    html_content <- create_interactive_html(enhanced_content, title)
    
    # Write final HTML
    writeLines(html_content, output_file)
    
    cat("✅ Interactive HTML generated successfully!\n")
    cat("🎯 Output saved to:", output_file, "\n")
    
    # Get file size
    file_size <- file.info(output_file)$size
    file_size_mb <- round(file_size / 1024 / 1024, 2)
    cat("📏 File size:", file_size_mb, "MB\n")
    
    cat("\n🎉 DOCUMENTATION GENERATION COMPLETE!\n")
    cat("Open", output_file, "in your web browser to view interactive flowcharts\n\n")
    
    return(TRUE)
    
  }, error = function(e) {
    cat("❌ ERROR generating HTML:", e$message, "\n")
    return(FALSE)
  })
}

#' Process Markdown Content for HTML
#'
#' Enhances markdown content for better HTML rendering
#' @param content Markdown content lines
#' @param title Document title
#' @return Enhanced markdown content
process_markdown_for_html <- function(content, title) {
  
  # Add flowchart numbering and anchors
  flowchart_counter <- 0
  enhanced_lines <- character()
  
  for (i in 1:length(content)) {
    line <- content[i]
    
    # Check if this is a flowchart header
    if (grepl("### \\*\\*Flowchart \\d+:", line)) {
      flowchart_counter <- flowchart_counter + 1
      
      # Extract flowchart title
      flowchart_title <- str_extract(line, "Flowchart \\d+: [^*]+")
      flowchart_title <- str_remove(flowchart_title, "\\*\\*$")
      
      # Create anchor ID
      anchor_id <- paste0("flowchart-", flowchart_counter)
      
      # Enhanced header with anchor and interactive features
      enhanced_header <- paste0(
        '<div class="flowchart-section" id="', anchor_id, '">',
        '\n### **', flowchart_title, '**',
        '\n<div class="flowchart-controls">',
        '\n  <button class="btn btn-sm btn-outline-primary" onclick="downloadFlowchart(\'', anchor_id, '\')">',
        '\n    📥 Download SVG</button>',
        '\n  <button class="btn btn-sm btn-outline-secondary" onclick="zoomFlowchart(\'', anchor_id, '\')">',
        '\n    🔍 Zoom</button>',
        '\n</div>'
      )
      
      enhanced_lines <- c(enhanced_lines, enhanced_header)
      
    } else if (grepl("```mermaid", line)) {
      # Add interactive wrapper for Mermaid diagrams
      diagram_id <- paste0("diagram-", flowchart_counter)
      enhanced_lines <- c(
        enhanced_lines,
        paste0('<div class="mermaid-wrapper" id="wrapper-', diagram_id, '">'),
        paste0('<div class="mermaid" id="', diagram_id, '">')
      )
      
    } else if (grepl("```$", line) && length(enhanced_lines) > 0 && grepl('<div class="mermaid"', tail(enhanced_lines, 1))) {
      # Close Mermaid wrapper
      enhanced_lines <- c(
        enhanced_lines,
        '</div>', # Close mermaid div
        '</div>', # Close wrapper div
        '</div>'  # Close flowchart section
      )
      
    } else {
      enhanced_lines <- c(enhanced_lines, line)
    }
  }
  
  return(enhanced_lines)
}

#' Create Interactive HTML
#'
#' Creates complete HTML document with Mermaid.js integration
#' @param content Processed markdown content
#' @param title Document title
#' @return Complete HTML content
create_interactive_html <- function(content, title) {
  
  # Convert markdown to HTML body
  temp_md <- tempfile(fileext = ".md")
  writeLines(content, temp_md)
  
  # Basic markdown to HTML conversion
  html_body <- markdown::markdownToHTML(temp_md, fragment.only = TRUE)
  
  # Clean up temp file
  unlink(temp_md)
  
  # Create complete HTML document
  html_template <- paste0(
    '<!DOCTYPE html>',
    '\n<html lang="en">',
    '\n<head>',
    '\n    <meta charset="UTF-8">',
    '\n    <meta name="viewport" content="width=device-width, initial-scale=1.0">',
    '\n    <title>', title, '</title>',
    '\n    ',
    '\n    <!-- Bootstrap CSS -->',
    '\n    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">',
    '\n    ',
    '\n    <!-- Mermaid.js -->',
    '\n    <script src="https://cdn.jsdelivr.net/npm/mermaid@10.6.1/dist/mermaid.min.js"></script>',
    '\n    ',
    '\n    <!-- Custom CSS -->',
    '\n    <style>',
    get_custom_css(),
    '\n    </style>',
    '\n</head>',
    '\n<body>',
    '\n    <!-- Navigation -->',
    '\n    <nav class="navbar navbar-expand-lg navbar-dark bg-primary fixed-top">',
    '\n        <div class="container-fluid">',
    '\n            <a class="navbar-brand" href="#top">',
    '\n                <strong>🧬 CAMK2D Pipeline Documentation</strong>',
    '\n            </a>',
    '\n            <div class="navbar-nav ms-auto">',
    '\n                <a class="btn btn-outline-light btn-sm me-2" href="CAMK_Analysis_Report.html" target="_blank">',
    '\n                    📈 Analysis Results',
    '\n                </a>',
    '\n                <button class="btn btn-outline-light btn-sm" onclick="toggleTOC()">',
    '\n                    📋 Table of Contents',
    '\n                </button>',
    '\n            </div>',
    '\n        </div>',
    '\n    </nav>',
    '\n    ',
    '\n    <!-- Table of Contents Sidebar -->',
    '\n    <div id="toc-sidebar" class="toc-sidebar">',
    '\n        <div class="toc-header">',
    '\n            <h5>📋 Contents</h5>',
    '\n            <button class="btn-close" onclick="toggleTOC()"></button>',
    '\n        </div>',
    '\n        <div class="toc-content">',
    generate_table_of_contents(content),
    '\n        </div>',
    '\n    </div>',
    '\n    ',
    '\n    <!-- Main Content -->',
    '\n    <div class="container-fluid main-content">',
    '\n        <div class="row">',
    '\n            <div class="col-12">',
    '\n                <div class="content-wrapper">',
    html_body,
    '\n                </div>',
    '\n            </div>',
    '\n        </div>',
    '\n    </div>',
    '\n    ',
    '\n    <!-- Bootstrap JS -->',
    '\n    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>',
    '\n    ',
    '\n    <!-- Custom JavaScript -->',
    '\n    <script>',
    get_custom_javascript(),
    '\n    </script>',
    '\n    ',
    '\n</body>',
    '\n</html>'
  )
  
  return(html_template)
}

#' Get Custom CSS
#'
#' Returns custom CSS for interactive documentation
#' @return CSS string
get_custom_css <- function() {
  css <- '
        body {
            padding-top: 70px;
            font-family: "Segoe UI", Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
        }
        
        .main-content {
            padding: 20px;
        }
        
        .content-wrapper {
            max-width: 1400px;
            margin: 0 auto;
        }
        
        /* Flowchart Sections */
        .flowchart-section {
            margin: 30px 0;
            padding: 20px;
            border: 1px solid #e0e0e0;
            border-radius: 8px;
            background: #fafafa;
        }
        
        .flowchart-controls {
            margin: 10px 0;
            text-align: right;
        }
        
        .flowchart-controls .btn {
            margin-left: 5px;
        }
        
        /* Mermaid Diagrams */
        .mermaid-wrapper {
            background: white;
            border: 1px solid #ddd;
            border-radius: 8px;
            padding: 20px;
            margin: 15px 0;
            overflow: auto;
        }
        
        .mermaid {
            text-align: center;
            font-family: "Segoe UI", sans-serif;
        }
        
        /* Table of Contents Sidebar */
        .toc-sidebar {
            position: fixed;
            top: 70px;
            left: -350px;
            width: 350px;
            height: calc(100vh - 70px);
            background: white;
            border-right: 1px solid #ddd;
            transition: left 0.3s ease;
            z-index: 1000;
            overflow-y: auto;
        }
        
        .toc-sidebar.open {
            left: 0;
        }
        
        .toc-header {
            padding: 15px 20px;
            border-bottom: 1px solid #eee;
            display: flex;
            justify-content: space-between;
            align-items: center;
            background: #f8f9fa;
        }
        
        .toc-content {
            padding: 20px;
        }
        
        .toc-item {
            display: block;
            padding: 8px 15px;
            text-decoration: none;
            color: #333;
            border-radius: 4px;
            margin: 2px 0;
            font-size: 0.9em;
            transition: background-color 0.2s;
        }
        
        .toc-item:hover {
            background-color: #e3f2fd;
            text-decoration: none;
        }
        
        .toc-section {
            font-weight: bold;
            color: #0d6efd;
            margin-top: 15px;
            padding: 5px 0;
            border-bottom: 1px solid #eee;
        }
        
        /* Responsive Design */
        @media (max-width: 768px) {
            .toc-sidebar {
                width: 100%;
                left: -100%;
            }
            
            .flowchart-controls {
                text-align: center;
            }
            
            .mermaid-wrapper {
                padding: 10px;
            }
        }
        
        /* Print Styles */
        @media print {
            .navbar, .toc-sidebar, .flowchart-controls {
                display: none !important;
            }
            
            body {
                padding-top: 0;
            }
            
            .flowchart-section {
                break-inside: avoid;
                margin: 20px 0;
            }
        }
        
        /* Code Blocks */
        pre {
            background: #f8f9fa;
            border: 1px solid #e9ecef;
            border-radius: 4px;
            padding: 15px;
            overflow-x: auto;
        }
        
        code {
            background: #f8f9fa;
            padding: 2px 4px;
            border-radius: 3px;
            font-size: 0.9em;
        }
        
        /* Tables */
        table {
            width: 100%;
            margin: 15px 0;
            border-collapse: collapse;
        }
        
        table th, table td {
            padding: 12px;
            text-align: left;
            border: 1px solid #ddd;
        }
        
        table th {
            background-color: #f8f9fa;
            font-weight: 600;
        }
        
        table tr:nth-child(even) {
            background-color: #f8f9fa;
        }
        
        /* Alert Boxes */
        .alert-custom {
            padding: 15px;
            margin: 15px 0;
            border-radius: 8px;
            border-left: 4px solid;
        }
        
        .alert-info {
            background-color: #cce7ff;
            border-left-color: #0066cc;
        }
        
        .alert-success {
            background-color: #d4edda;
            border-left-color: #28a745;
        }
        
        .alert-warning {
            background-color: #fff3cd;
            border-left-color: #ffc107;
        }
  '
  
  return(css)
}

#' Generate Table of Contents
#'
#' Creates navigation menu for all flowcharts
#' @param content Markdown content
#' @return HTML for table of contents
generate_table_of_contents <- function(content) {
  
  toc_html <- character()
  flowchart_counter <- 0
  current_section <- ""
  
  for (line in content) {
    # Check for main sections
    if (grepl("^## ", line)) {
      section_title <- str_remove(line, "^## \\*\\*|\\*\\*$")
      section_title <- str_remove_all(section_title, "[🏗️📊🔬🛡️⚡📋]")
      section_title <- str_trim(section_title)
      current_section <- section_title
      
      toc_html <- c(toc_html, paste0('<div class="toc-section">', section_title, '</div>'))
    }
    
    # Check for flowcharts
    if (grepl("### \\*\\*Flowchart \\d+:", line)) {
      flowchart_counter <- flowchart_counter + 1
      
      # Extract flowchart title
      flowchart_title <- str_extract(line, "Flowchart \\d+: [^*]+")
      flowchart_title <- str_remove(flowchart_title, "\\*\\*$")
      flowchart_title <- str_trim(flowchart_title)
      
      # Create TOC link
      anchor_id <- paste0("flowchart-", flowchart_counter)
      toc_link <- paste0(
        '<a href="#', anchor_id, '" class="toc-item" onclick="closeTOC()">',
        '📊 ', flowchart_title,
        '</a>'
      )
      
      toc_html <- c(toc_html, toc_link)
    }
  }
  
  return(paste(toc_html, collapse = '\n'))
}

#' Get Custom JavaScript
#'
#' Returns custom JavaScript for interactivity
#' @return JavaScript string
get_custom_javascript <- function() {
  js <- '
        // Initialize Mermaid
        mermaid.initialize({
            startOnLoad: true,
            theme: "default",
            themeVariables: {
                primaryColor: "#e3f2fd",
                primaryTextColor: "#1976d2",
                primaryBorderColor: "#1976d2",
                lineColor: "#666",
                secondaryColor: "#f5f5f5",
                tertiaryColor: "#ffffff"
            },
            flowchart: {
                useMaxWidth: true,
                htmlLabels: true,
                curve: "cardinal"
            },
            securityLevel: "loose",
            maxTextSize: 90000,
            maxEdges: 2000,
            deterministicIds: true,
            deterministicIDSeed: "mermaidFlowchart"
        });
        
        // Error handling for Mermaid rendering
        mermaid.parseError = function(err, hash) {
            console.warn("Mermaid parsing error:", err);
            if (hash && hash.str) {
                const errorDiv = document.createElement("div");
                errorDiv.className = "alert alert-warning";
                errorDiv.innerHTML = "⚠️ Flowchart too complex to render. Please simplify or split into smaller diagrams.";
                hash.str.parentNode.replaceChild(errorDiv, hash.str);
            }
        };
        
        // Table of Contents Toggle
        function toggleTOC() {
            const sidebar = document.getElementById("toc-sidebar");
            sidebar.classList.toggle("open");
        }
        
        function closeTOC() {
            const sidebar = document.getElementById("toc-sidebar");
            sidebar.classList.remove("open");
        }
        
        // Download Flowchart as SVG
        function downloadFlowchart(flowchartId) {
            const diagramId = flowchartId.replace("flowchart-", "diagram-");
            const svgElement = document.querySelector(`#${diagramId} svg`);
            
            if (svgElement) {
                const svgData = new XMLSerializer().serializeToString(svgElement);
                const svgBlob = new Blob([svgData], {type: "image/svg+xml;charset=utf-8"});
                const svgUrl = URL.createObjectURL(svgBlob);
                
                const downloadLink = document.createElement("a");
                downloadLink.href = svgUrl;
                downloadLink.download = `${flowchartId}.svg`;
                document.body.appendChild(downloadLink);
                downloadLink.click();
                document.body.removeChild(downloadLink);
                URL.revokeObjectURL(svgUrl);
                
                // Success notification
                showNotification("Flowchart downloaded successfully!", "success");
            } else {
                showNotification("Unable to find flowchart for download", "error");
            }
        }
        
        // Zoom Flowchart
        function zoomFlowchart(flowchartId) {
            const diagramId = flowchartId.replace("flowchart-", "diagram-");
            const wrapper = document.getElementById(`wrapper-${diagramId}`);
            
            if (wrapper) {
                wrapper.classList.toggle("zoomed");
                
                if (wrapper.classList.contains("zoomed")) {
                    wrapper.style.transform = "scale(1.5)";
                    wrapper.style.transformOrigin = "top center";
                    wrapper.style.margin = "50px 0";
                } else {
                    wrapper.style.transform = "scale(1)";
                    wrapper.style.margin = "15px 0";
                }
            }
        }
        
        // Show Notification
        function showNotification(message, type = "info") {
            const notification = document.createElement("div");
            notification.className = `alert alert-${type} position-fixed`;
            notification.style.top = "80px";
            notification.style.right = "20px";
            notification.style.zIndex = "9999";
            notification.style.minWidth = "300px";
            notification.innerHTML = message;
            
            document.body.appendChild(notification);
            
            setTimeout(() => {
                notification.remove();
            }, 3000);
        }
        
        // Smooth Scrolling for TOC Links
        document.addEventListener("DOMContentLoaded", function() {
            const tocLinks = document.querySelectorAll(".toc-item");
            
            tocLinks.forEach(link => {
                link.addEventListener("click", function(e) {
                    e.preventDefault();
                    const targetId = this.getAttribute("href").substring(1);
                    const targetElement = document.getElementById(targetId);
                    
                    if (targetElement) {
                        targetElement.scrollIntoView({
                            behavior: "smooth",
                            block: "start"
                        });
                    }
                });
            });
        });
        
        // Keyboard Shortcuts
        document.addEventListener("keydown", function(e) {
            // Ctrl/Cmd + T to toggle TOC
            if ((e.ctrlKey || e.metaKey) && e.key === "t") {
                e.preventDefault();
                toggleTOC();
            }
            
            // Escape to close TOC
            if (e.key === "Escape") {
                closeTOC();
            }
        });
        
        // Auto-generate anchor links for sections
        document.addEventListener("DOMContentLoaded", function() {
            const headings = document.querySelectorAll("h1, h2, h3, h4");
            
            headings.forEach((heading, index) => {
                if (!heading.id) {
                    const text = heading.textContent.toLowerCase()
                        .replace(/[^\\w\\s-]/g, "")
                        .replace(/\\s+/g, "-");
                    heading.id = text || `heading-${index}`;
                }
            });
        });
        
        // Print Functionality
        function printDocumentation() {
            window.print();
        }
        
        // Search Functionality (basic)
        function searchFlowcharts() {
            const searchTerm = prompt("Search for flowchart:");
            if (searchTerm) {
                const tocItems = document.querySelectorAll(".toc-item");
                let found = false;
                
                tocItems.forEach(item => {
                    if (item.textContent.toLowerCase().includes(searchTerm.toLowerCase())) {
                        item.click();
                        found = true;
                        return;
                    }
                });
                
                if (!found) {
                    showNotification(`No flowchart found containing "${searchTerm}"`, "warning");
                }
            }
        }
  '
  
  return(js)
}

#' Generate Analysis Results Summary
#'
#' Creates a summary section with current analysis results
#' @return Character vector with markdown content
generate_results_summary <- function() {
  summary_content <- character()
  
  # Try to load analysis results
  meta_file <- "output/current/CAMK_meta_analysis_FINAL.csv"
  dge_file <- "output/current/CAMK_DGE_all_6_datasets_COMPREHENSIVE.csv"
  
  summary_content <- c(summary_content, 
    "## 📊 **Current Analysis Results**",
    "",
    '<div class="alert-custom alert-info">',
    "**Latest Pipeline Execution Results** - Live data from current analysis run",
    "</div>",
    ""
  )
  
  if (file.exists(meta_file)) {
    tryCatch({
      meta_results <- read.csv(meta_file, stringsAsFactors = FALSE)
      
      # Find CAMK2D results
      camk2d_row <- meta_results[meta_results$Gene == "CAMK2D", ]
      
      if (nrow(camk2d_row) > 0) {
        summary_content <- c(summary_content,
          "### **🎯 Primary Gene Results**",
          "",
          paste("- **CAMK2D Log Fold Change:** ", round(camk2d_row$logFC, 4)),
          paste("- **P-value:** ", format(camk2d_row$P.Value, scientific = TRUE, digits = 3)),
          paste("- **Significance:** ", ifelse(camk2d_row$P.Value < 0.05, "✅ Significant", "❌ Not Significant")),
          "",
          "### **📈 Meta-Analysis Summary**",
          "",
          paste("- **Total Genes Analyzed:** ", nrow(meta_results)),
          paste("- **Significant Genes:** ", sum(meta_results$P.Value < 0.05, na.rm = TRUE)),
          paste("- **Effect Direction:** ", ifelse(camk2d_row$logFC > 0, "Upregulated", "Downregulated")),
          ""
        )
      }
    }, error = function(e) {
      summary_content <- c(summary_content, "⚠️ Meta-analysis results not yet available", "")
    })
  } else {
    summary_content <- c(summary_content, "ℹ️ Meta-analysis results will appear here after pipeline execution", "")
  }
  
  # Add link to full analysis report
  summary_content <- c(summary_content,
    "### **🔗 Quick Links**",
    "",
    "- [📊 View Complete Analysis Report](CAMK_Analysis_Report.html)",
    "- [📈 Download Meta-Analysis Results](CAMK_meta_analysis_FINAL.csv)",
    "- [📋 Download DGE Results](CAMK_DGE_all_6_datasets_COMPREHENSIVE.csv)",
    "",
    "---",
    ""
  )
  
  return(summary_content)
}

# Main execution
if (!interactive()) {
  cat("🚀 Starting Interactive Documentation Generator...\n")
  
  result <- generate_interactive_documentation()
  
  if (result) {
    cat("🎉 SUCCESS: Interactive documentation generated!\n")
    cat("📱 Open 'Interactive_Technical_Documentation.html' in your browser\n")
  } else {
    cat("❌ FAILED: Documentation generation failed\n")
    quit(status = 1)
  }
}

cat("✅ Interactive Documentation Generator loaded successfully\n\n")