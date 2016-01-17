function kksendmail_janelia(mail_addr, subject, msg, attachment)

setpref('Internet','SMTP_Server','mail.janelia.org');
setpref('Internet','E_mail','khairyk@janelia.hhmi.org');

if nargin == 0,sendmail('khairyk@janelia.hhmi.org','***MATLAB***');
elseif nargin==2
    subject = ['***MATLAB*** ' subject];
  sendmail(mail_addr,subject);
elseif nargin==3
    subject = ['***MATLAB*** ' subject];
  sendmail(mail_addr,subject, msg);
elseif nargin==4
    subject = ['***MATLAB*** ' subject];
  sendmail(mail_addr,subject, msg, {attachment});
end

