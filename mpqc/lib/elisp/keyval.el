
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
  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults '(keyval-font-lock-keywords t t))
  (run-hooks 'keyval-mode-hook))

(defun keyval-setup ())

(defvar keyval-key-face 'keyval-key-face
  "Face for keys in KeyVal file.")
(make-face keyval-key-face)
(make-face-bold keyval-key-face)
(set-face-foreground keyval-key-face "Cyan")

(defvar keyval-classtype-face 'keyval-classtype-face
  "Face for names of classes in KeyVal file.")
(make-face keyval-classtype-face)
(set-face-foreground keyval-classtype-face "Green")

(defvar keyval-reference-face 'keyval-reference-face
  "Face for references in KeyVal file.")
(make-face keyval-reference-face)
(set-face-foreground keyval-reference-face "Orange")
(set-face-underline-p keyval-reference-face t)

(defvar keyval-font-lock-keywords
  '(("%.*" . font-lock-comment-face)
    ("<.*>" . keyval-classtype-face)
    ("\"[^\"\n]+\"" . font-lock-string-face)
    ("$[A-Za-z0-9_\.:*+-/]*" . keyval-reference-face)
    ("{ *\\([A-Za-z0-9_\.*+-/ ]*\\>\\) *} *=" (1 keyval-key-face))
    ("\\([A-Za-z0-9_\.*+-/]*\\>\\) *=" (1 keyval-key-face))
    )
  "Default expressions to highlight in KeyVal mode.")
