
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
  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults '(mpqc-font-lock-keywords t))
  (run-hooks 'mpqc-mode-hook))

(defun mpqc-setup ()
  (toggle-read-only 1)
  )

(defvar mpqc-success-face 'mpqc-success-face
  "Face for usable mpqc output.")
(make-face mpqc-success-face)
(make-face-bold mpqc-success-face)
(set-face-foreground mpqc-success-face "Cyan")

(defvar mpqc-info-face 'mpqc-info-face
  "Face for informational mpqc output.")
(make-face mpqc-info-face)
(make-face-bold mpqc-info-face)
(set-face-foreground mpqc-info-face "Green")

(defvar mpqc-warning-face 'mpqc-warning-face
  "Face for warnings in MPQC output.")
(make-face mpqc-warning-face)
(set-face-foreground mpqc-warning-face "Orange")
(set-face-underline-p mpqc-warning-face t)

(defvar mpqc-font-lock-keywords
  '(("converged scf energy.*" . mpqc-success-face)
    ("Converged.*Internal Coordinates" . mpqc-success-face)
    ("Nonconverged.*Internal Coordinates" . mpqc-warning-face)
    ("Fixed.*Internal Coordinates" . mpqc-warning-face)
    ("Initial.*Internal Coordinates" . mpqc-info-face)
    ("Updated.*Internal Coordinates" . mpqc-info-face)
    ("internal coordinates" . mpqc-info-face)
    ("max of 1/2 idisp.*" . mpqc-info-face)
    ("DEBUG.*" . mpqc-warning-face)
    ("WARNING.*" . mpqc-warning-face)
    ("NOTICE.*" . mpqc-warning-face)
    ("Too many geometry.*" . mpqc-warning-face)
    ("TORS" . mpqc-info-face)
    ("BEND" . mpqc-info-face)
    ("STRE" . mpqc-info-face)
    ("OUTP" . mpqc-info-face)
    )
  "Default expressions to highlight in MPQC mode.")
