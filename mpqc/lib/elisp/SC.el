;;
;; stuff for formatting files in the SC libraries
;;

(require 'cc-mode)

(setq c-basic-half-offset 2)

(defun clj-adaptive-block-open (langelem)
  ;; when substatement is on semantics list, return
  ;; -(c-basic-offset - c-basic-half-offset) to give a
  ;; total offset of c-basic-half-offset,
  ;; otherwise return c-basic-half-offset
  (if (assq 'substatement c-semantics)
      (+ c-basic-half-offset (- c-basic-offset))
    c-basic-half-offset))

(defun clj-lineup-math (langelem)
  ;; line up math statement-cont so that stuff after the "+", "-", etc
  ;; lines up with the stuff after the equals
  (save-excursion
    (let ((adjustment (progn
                        (beginning-of-line)
                        (skip-chars-forward " \t" (c-point 'eol))
                        (- (current-column)
                           (progn (skip-chars-forward " \t+-/*" (c-point 'eol))
                                  (current-column)))))
          (curcol (progn
		    (goto-char (cdr langelem))
		    (current-column))))
      (skip-chars-forward "^=" (c-point 'eol))
      (if (/= (following-char) ?=)
	  ;; there's no equal sign on the line
	  c-basic-offset
	;; calculate indentation column after equals and ws and sign
	(forward-char 1)
	(skip-chars-forward " \t-")
	(+ (- (current-column) curcol) adjustment))
      )))

(defun clj-adaptive-block-close (langelem)
  ;; these closes blocks in a way that is consistent with the way
  ;; clj-adaptive-statement-block-intro indents the first statement
  (- (clj-adaptive-statement-block-intro langelem)
     (- c-basic-offset c-basic-half-offset))
)

(defun clj-adaptive-statement-block-intro (langelem)
  ;; this lines up the first statement in a block by a full basic
  ;; offset, unless we are lining up to a "{" which is already
  ;; half indented
  (save-excursion
    (progn
      (goto-char (cdr langelem))
      (if (/= (following-char) ?{)
          ;; next char is not a "{"
          c-basic-offset
        ;; use remainder of half offset
        (- c-basic-offset c-basic-half-offset))
      )))

;;
;; this is the style to use when editting Ed's files
;;
(c-add-style "ETS" '((c-basic-offset . 2)
                     (c-offsets-alist . ((access-label      . -)
                                         (inclass           . ++)
                                         (label             . 0)
                                         ))
                     )
             )

;;
;; this is the style to use when editting Curt's files
;;
(c-add-style "CLJ" '(
    (c-offsets-alist . (
        (block-open      . clj-adaptive-block-open)
        (statement       . c-lineup-runin-statements)
        (statement-cont  . clj-lineup-math)
        (statement-block-intro . clj-adaptive-statement-block-intro)
        (defun-block-intro . 2)
        (inher-intro . 2)
        (access-label . -2)
        (block-close . clj-adaptive-block-close)
        (member-init-intro . 2)
        )
    ))
)

(defun clj-style ()
  "Change to insane C indentation"
  (interactive)
  (set-c-style "CLJ")
  )
(defun ets-style ()
  "Change to sensible C indentation"
  (interactive)
  (set-c-style "ETS")
  )

(define-key c-mode-map "\C-ce" 'ets-style)
(define-key c-mode-map "\C-cj" 'clj-style)
(define-key c-mode-map "\C-j"  'reindent-then-newline-and-indent)
(define-key c-mode-map "\C-m"  'newline-and-indent)

;;
;; stuff for CLJ's compile hacks
;;

(defun compile-modify-path (thisdir)
  (let ((tmpdir (expand-file-name thisdir)))
    (setq thisdir "")
    (while (>= (length tmpdir) (length sc-src-dir))
      (if (string= (substring tmpdir 0 (length sc-src-dir)) sc-src-dir)
          (let ()
            (setq thisdir (concat thisdir sc-arch-dir))
            (setq tmpdir (substring tmpdir (length sc-src-dir) nil))
            )
        (let ()
          (setq thisdir (concat thisdir (substring tmpdir 0 1)))
          (setq tmpdir (substring tmpdir 1 nil))
          )
        )
      )
    (setq thisdir (concat thisdir tmpdir))
    )
  thisdir
)
