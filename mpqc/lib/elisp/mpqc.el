
;;;###autoload
(defvar mpqc-mode-hook nil
  "*List of hook functions run by `mpqc-mode' (see `run-hooks').")

(defvar mpqc-mode-map
  (let ((map (make-sparse-keymap)))
    (define-key map " " 'scroll-up)
    (define-key map "\^?" 'scroll-down)
    map)
  "Keymap for mpqc output buffers.")

(defun mpqc-mode ()
  "Major mode for mpqc output files."
  (interactive)
  (fundamental-mode)
  (use-local-map mpqc-mode-map)
  (setq major-mode 'mpqc-mode
	mode-name "Mpqc")
  (mpqc-setup)
  (run-hooks 'mpqc-mode-hook))

(defun mpqc-setup ()
  (toggle-read-only 1)
  )
