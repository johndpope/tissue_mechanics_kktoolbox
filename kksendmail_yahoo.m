function kksendmail_yahoo(mail_addr, subject, msg, attachment)
%% example:
%%% kksendmail_yahoo('khaledkhairy@yahoo.com', 'test','testing msg','ds.m');
myaddress = 'khaledkhairy@yahoo.com';
mypassword = 'apache0455';

setpref('Internet','E_mail',myaddress);
setpref('Internet','SMTP_Server','smtp.mail.yahoo.com');
setpref('Internet','SMTP_Username',myaddress);
setpref('Internet','SMTP_Password',mypassword);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
                  'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');


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

