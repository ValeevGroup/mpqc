
;;;###autoload
(defvar keyval-mode-hook nil
  "*List of hook functions run by `keyval-mode' (see `run-hooks').")

(defvar keyval-mode-map
  (let ((map (make-sparse-keymap)))
    ;(define-key map " " 'scroll-up)
    ;(define-key map "\^?" 'scroll-down)
    map)
  "Keymap for KeyVal input buffers.")

(defun keyval-mode ()
  "Major mode for KeyVal input files."
  (interactive)
  (fundamental-mode)
  (use-local-map keyval-mode-map)
  (setq major-mode 'keyval-mode
	mode-name "KeyVal")
  (keyval-setup)
  (run-hooks 'keyval-mode-hook))

(defun keyval-setup ()
  )
