use std::fmt::{Debug, Display, Formatter};
use std::num::ParseFloatError;
use std::sync::mpsc::{RecvError, RecvTimeoutError, SendError};
use std::time::SystemTimeError;

mod names {
    pub const MOCASA: &str = "Mocasa error";
    pub const IO: &str = "I/O error";
    pub const TOML_DE: &str = "TOML deserialization error";
    pub const TOML_SER: &str = "TOML serialization error";
    pub const PARSE_FLOAT: &str = "parse float error";
    pub const SEND: &str = "send error";
    pub const RECEIVE: &str = "receive error";
    pub const RECEIVE_TIMEOUT: &str = "receive timeout error";
    pub const SYSTEM_TIME: &str = "system time error";
    pub const SERDE_JSON: &str = "Serde JSON";
}

pub enum ErrorKind {
    Mocasa,
    IOError,
    TomlDe,
    TomlSer,
    ParseFloat,
    Send,
    Receive,
    ReceiveTimeout,
    SystemTime,
    SerdeJson,
}

pub struct Error {
    kind: ErrorKind,
    message: String,
}

impl From<&str> for Error {
    fn from(message: &str) -> Self {
        Error::new(ErrorKind::Mocasa, message.to_string())
    }
}

impl From<String> for Error {
    fn from(message: String) -> Self {
        Error::new(ErrorKind::Mocasa, message)
    }
}

impl From<std::io::Error> for Error {
    fn from(io_error: std::io::Error) -> Self {
        let message = io_error.to_string();
        Error::new(ErrorKind::IOError, message)
    }
}

impl From<toml::de::Error> for Error {
    fn from(toml_de_error: toml::de::Error) -> Self {
        let message = toml_de_error.to_string();
        Error::new(ErrorKind::TomlDe, message)
    }
}

impl From<toml::ser::Error> for Error {
    fn from(toml_ser_error: toml::ser::Error) -> Self {
        let message = toml_ser_error.to_string();
        Error::new(ErrorKind::TomlSer, message)
    }
}

impl From<ParseFloatError> for Error {
    fn from(parse_float_error: ParseFloatError) -> Self {
        let message = parse_float_error.to_string();
        Error::new(ErrorKind::ParseFloat, message)
    }
}

impl<O> From<SendError<O>> for Error {
    fn from(send_error: SendError<O>) -> Self {
        let message = send_error.to_string();
        Error::new(ErrorKind::Send, message)
    }
}

impl From<RecvError> for Error {
    fn from(receive_error: RecvError) -> Self {
        let message = receive_error.to_string();
        Error::new(ErrorKind::Receive, message)
    }
}

impl From<RecvTimeoutError> for Error {
    fn from(receive_timeout_error: RecvTimeoutError) -> Self {
        let message = receive_timeout_error.to_string();
        Error::new(ErrorKind::ReceiveTimeout, message)
    }
}

impl From<SystemTimeError> for Error {
    fn from(system_time_error: SystemTimeError) -> Self {
        let message = system_time_error.to_string();
        Error::new(ErrorKind::SystemTime, message)
    }
}

impl From<serde_json::Error> for Error {
    fn from(serde_json_error: serde_json::Error) -> Self {
        let message = serde_json_error.to_string();
        Error::new(ErrorKind::SerdeJson, message)
    }
}

impl ErrorKind {
    pub fn as_str(&self) -> &'static str {
        match self {
            ErrorKind::Mocasa => { names::MOCASA }
            ErrorKind::IOError => { names::IO }
            ErrorKind::TomlDe => { names::TOML_DE }
            ErrorKind::TomlSer => { names::TOML_SER }
            ErrorKind::ParseFloat => { names::PARSE_FLOAT }
            ErrorKind::Send => { names::SEND }
            ErrorKind::Receive => { names::RECEIVE }
            ErrorKind::ReceiveTimeout => { names::RECEIVE_TIMEOUT }
            ErrorKind::SystemTime => { names::SYSTEM_TIME }
            ErrorKind::SerdeJson => { names::SERDE_JSON }
        }
    }
}

impl Display for ErrorKind {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

impl Error {
    pub(crate) fn new(kind: ErrorKind, message: String) -> Error {
        Error { kind, message }
    }
}

impl Debug for Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.kind.as_str(), self.message)
    }
}

impl Display for Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.kind.as_str(), self.message)
    }
}

impl std::error::Error for Error {}

pub(crate) fn for_file<T>(file: &str, result: std::io::Result<T>) -> Result<T, Error> {
    match result {
        Ok(value) => { Ok(value) }
        Err(error) => {
            Err(Error::new(ErrorKind::IOError, format!("{}: {}", file, error)))
        }
    }
}