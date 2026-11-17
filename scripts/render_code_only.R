# Create a code-and-explanations page from the step-by-step R Markdown.
lines <- readLines("workflow.Rmd", warn = FALSE, encoding = "UTF-8")
escape <- function(x) {
 x <- gsub("&", "&amp;", x, fixed = TRUE)
 x <- gsub("<", "&lt;", x, fixed = TRUE)
 gsub(">", "&gt;", x, fixed = TRUE)
}
parts <- c('<!doctype html><html lang="en"><head><meta charset="utf-8">',
 '<meta name="viewport" content="width=device-width, initial-scale=1"><title>Workshop commands</title>',
 '<style>body{font-family:system-ui,sans-serif;margin:2rem auto;max-width:1000px;padding:0 1rem;color:#223344;line-height:1.5}',
 'h1,h2{color:#14517a}h2{margin-top:2rem;border-bottom:1px solid #d8e4ea}pre{overflow:auto;background:#f2f6f8;border-left:4px solid #4a91aa;border-radius:4px;padding:2.8rem 1rem 1rem;line-height:1.45}',
 'code{font-family:ui-monospace,monospace}.code-block{position:relative}.copy-button{position:absolute;right:.6rem;top:.55rem;padding:.3rem .65rem;border:1px solid #7495a4;border-radius:4px;background:#fff;color:#17465d;cursor:pointer}',
 '.copy-button:focus-visible{outline:3px solid #217da8;outline-offset:2px}</style></head><body>',
 '<h1>Workshop commands</h1><p>Open the workshop folder in RStudio. Work from top to bottom; click Copy and paste each complete command into the Console. This page contains no analysis results.</p>')
starts <- grep('^```\\{r step-', lines)
for (start in starts) {
  previous_headings <- which(grepl('^#{2,3} ', lines[seq_len(start - 1L)]))
  if (length(previous_headings) && (!exists('last_heading') || tail(previous_headings, 1L) != last_heading)) {
    last_heading <- tail(previous_headings, 1L)
    level <- if (startsWith(lines[last_heading], '### ')) 'h3' else 'h2'
    parts <- c(parts, paste0('<', level, '>', escape(sub('^#{2,3} ', '', lines[last_heading])), '</', level, '>'))
  }
  end <- which(seq_along(lines) > start & lines == '```')[1]
  description <- lines[start - 2L]
  code <- lines[seq.int(start + 1L, end - 1L)]
  parts <- c(parts, paste0('<p>', escape(description), '</p><div class="code-block"><button class="copy-button" type="button" aria-label="Copy this code block">Copy</button><pre><code>'), escape(code), '</code></pre></div>')
}
parts <- c(parts, '<script>',
 'document.addEventListener("click", async function(event) {',
 'const button=event.target.closest(".copy-button"); if(!button)return;',
 'const code=button.parentElement.querySelector("code").textContent.replace(/^\\n|\\n$/g, ""); let success=false;',
 'try { if(navigator.clipboard && window.isSecureContext){await navigator.clipboard.writeText(code);success=true;} } catch(error){}',
 'if(!success){const helper=document.createElement("textarea");helper.value=code;helper.setAttribute("readonly", "");helper.style.position="fixed";helper.style.opacity="0";document.body.appendChild(helper);helper.select();try{success=document.execCommand("copy");}catch(error){}helper.remove();}',
 'button.textContent=success?"Copied!":"Select code to copy";setTimeout(function(){button.textContent="Copy";},2500);',
 '});</script></body></html>')
writeLines(parts, "workflow-code.html", useBytes = TRUE)
