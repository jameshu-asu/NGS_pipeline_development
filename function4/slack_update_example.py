from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError

def slack_bot(channel:str, message:str) -> None:
    SLACK_BOT_TOKEN = 'xoxb-slack-token'
    client = WebClient(token=SLACK_BOT_TOKEN)
    result = client.chat_postMessage(channel=channel, text=message)
    return print(result)

slack_bot('C04PVD6D1ST', 'Hello World!')









